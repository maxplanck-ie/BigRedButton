import configparser
import json
from unittest.mock import Mock, patch

import requests

from BRB.misc import resolveDeliverTo, resolveGroup


def _make_config(deliver_to=None):
    config = configparser.ConfigParser()
    config["Parkour"] = {
        "InternalPIsURL": "https://parkour-demo.example.org/api/internal_pis/",
        "user": "jefke",
        "password": "123",
        "cert": "",
    }
    config["Internals"] = {"Organizations": "MPI-IE"}
    if deliver_to is not None:
        config["Internals"]["deliverTo"] = json.dumps(deliver_to)
    return config


class TestResolveDeliverTo:
    @patch("BRB.misc.requests.get")
    def test_returns_only_non_empty_overrides(self, mock_get):
        mock_resp = Mock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {
            "pis": {
                "cabezas-wallscheid": "cabezas",
                "akhtar": None,
                "alhaj abed": "alhajabed",
            }
        }
        mock_get.return_value = mock_resp
        config = _make_config()

        result = resolveDeliverTo(config)

        assert result == {"cabezas-wallscheid": "cabezas", "alhaj abed": "alhajabed"}
        mock_get.assert_called_once_with(
            "https://parkour-demo.example.org/api/internal_pis/",
            params={"organizations": "MPI-IE"},
            auth=("jefke", "123"),
            verify="",
        )

    @patch("BRB.misc.requests.get")
    def test_network_failure_falls_back_to_empty_map(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("boom")
        config = _make_config()

        assert resolveDeliverTo(config) == {}

    @patch("BRB.misc.requests.get")
    def test_non_200_falls_back_to_empty_map(self, mock_get):
        mock_resp = Mock()
        mock_resp.status_code = 500
        mock_resp.raise_for_status.side_effect = requests.exceptions.HTTPError(
            "500 error"
        )
        mock_get.return_value = mock_resp
        config = _make_config()

        assert resolveDeliverTo(config) == {}

    @patch("BRB.misc.requests.get")
    def test_legacy_bare_list_response_has_no_overrides(self, mock_get):
        mock_resp = Mock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {"pis": ["akhtar", "cabezas-wallscheid"]}
        mock_get.return_value = mock_resp
        config = _make_config()

        assert resolveDeliverTo(config) == {}


class TestResolveGroup:
    def test_uses_deliver_to_override_when_present(self):
        config = _make_config(deliver_to={"cabezas-wallscheid": "cabezas"})
        assert resolveGroup(config, "4035_Demollin_Cabezas-Wallscheid") == "cabezas"

    def test_falls_back_to_hyphen_truncation_without_override(self):
        config = _make_config(deliver_to={})
        assert resolveGroup(config, "4035_Demollin_Cabezas-Wallscheid") == "cabezas"

    def test_falls_back_when_deliverTo_key_missing_entirely(self):
        config = _make_config()  # no "Internals"/"deliverTo" key set at all
        assert resolveGroup(config, "4035_Demollin_Cabezas-Wallscheid") == "cabezas"

    def test_simple_surname_unaffected(self):
        config = _make_config(deliver_to={})
        assert resolveGroup(config, "1234_jdoe_manke") == "manke"

    def test_override_wins_even_when_it_differs_from_fallback_guess(self):
        # A PI whose deliver_to override isn't just "truncate at the first
        # hyphen" - the override must take priority over the guess.
        config = _make_config(deliver_to={"alhaj abed": "alhajabed"})
        assert resolveGroup(config, "9001_jdoe_AlHaj Abed") == "alhajabed"
