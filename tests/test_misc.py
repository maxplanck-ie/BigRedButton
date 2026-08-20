import configparser
import json
import subprocess as sp
from unittest.mock import Mock, mock_open, patch

import pytest
import requests

from BRB import misc
from BRB.misc import resolveDeliverTo, resolveGroup


class TestLoadUserDictionary:
    def test_parses_tab_separated_user_file(self):
        content = "loginA\tnameA\temailA\nloginB\tnameB\temailB\n"
        with patch("builtins.open", mock_open(read_data=content)):
            d = misc.loadUserDictionary()
        assert d == {
            "nameA": ["loginA", "emailA\n"],
            "nameB": ["loginB", "emailB\n"],
        }


class TestGetLatestSeqdir:
    def test_no_sequencing_data_dirs_returns_bare_name(self, tmp_path, monkeypatch):
        piDir = tmp_path / "PI"
        piDir.mkdir()
        (piDir / "other_stuff").mkdir()
        assert misc.getLatestSeqdir(str(tmp_path), "PI") == "sequencing_data"

    def test_bare_sequencing_data_dir_returns_bare_name(self, tmp_path):
        piDir = tmp_path / "PI"
        piDir.mkdir()
        (piDir / "sequencing_data").mkdir()
        assert misc.getLatestSeqdir(str(tmp_path), "PI") == "sequencing_data"

    def test_picks_highest_numbered_suffix(self, tmp_path):
        piDir = tmp_path / "PI"
        piDir.mkdir()
        (piDir / "sequencing_data").mkdir()
        (piDir / "sequencing_data2").mkdir()
        (piDir / "sequencing_data10").mkdir()
        assert misc.getLatestSeqdir(str(tmp_path), "PI") == "sequencing_data10"


class TestFormatGitDescribe:
    def test_clean_tag_is_reformatted(self):
        assert misc._formatGitDescribe("v1.2.3-4-gabc1234") == "v1.2.3 +4:gabc1234"

    def test_dirty_tag_keeps_dirty_suffix(self):
        assert (
            misc._formatGitDescribe("v1.2.3-4-gabc1234-dirty")
            == "v1.2.3 +4:gabc1234-dirty"
        )

    def test_bare_hash_with_no_tag_is_returned_unchanged(self):
        assert misc._formatGitDescribe("abc1234") == "abc1234"


class TestGetVersion:
    def test_uses_git_describe_when_available(self, monkeypatch):
        def fakeRun(cmd, **kwargs):
            return sp.CompletedProcess(cmd, 0, stdout="v1.0.0-2-gdeadbee\n")

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        assert misc.getVersion("BRB") == "v1.0.0 +2:gdeadbee"

    def test_falls_back_to_installed_metadata_when_git_missing(self, monkeypatch):
        def fakeRun(cmd, **kwargs):
            raise FileNotFoundError

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        monkeypatch.setattr(misc, "version", lambda distName: "9.9.9")
        assert misc.getVersion("BRB") == "9.9.9"

    def test_falls_back_when_not_a_git_repo(self, monkeypatch):
        def fakeRun(cmd, **kwargs):
            raise sp.CalledProcessError(128, cmd)

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        monkeypatch.setattr(misc, "version", lambda distName: "9.9.9")
        assert misc.getVersion("BRB") == "9.9.9"

    def test_falls_back_on_timeout(self, monkeypatch):
        def fakeRun(cmd, **kwargs):
            raise sp.TimeoutExpired(cmd, 5)

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        monkeypatch.setattr(misc, "version", lambda distName: "9.9.9")
        assert misc.getVersion("BRB") == "9.9.9"


class TestConfigGitInfo:
    def test_returns_none_when_not_inside_a_git_repo(self, tmp_path, monkeypatch):
        def fakeRun(cmd, **kwargs):
            raise sp.CalledProcessError(128, cmd)

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        assert misc.configGitInfo(str(tmp_path / "brb.ini")) is None

    def test_returns_none_when_git_binary_missing(self, tmp_path, monkeypatch):
        def fakeRun(cmd, **kwargs):
            raise FileNotFoundError

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        assert misc.configGitInfo(str(tmp_path / "brb.ini")) is None

    def test_clean_tree_returns_commit_hash(self, tmp_path, monkeypatch):
        calls = []

        def fakeRun(cmd, **kwargs):
            calls.append(cmd)
            if "rev-parse" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="true\n")
            if "status" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="")
            if "log" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="abc1234\n")
            raise AssertionError(f"unexpected command {cmd}")

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        result = misc.configGitInfo(str(tmp_path / "brb.ini"))
        assert result == "abc1234"

    def test_dirty_tree_exits(self, tmp_path, monkeypatch, capsys):
        def fakeRun(cmd, **kwargs):
            if "rev-parse" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="true\n")
            if "status" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout=" M brb.ini\n")
            raise AssertionError(f"unexpected command {cmd}")

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        with pytest.raises(SystemExit) as excinfo:
            misc.configGitInfo(str(tmp_path / "brb.ini"))
        assert excinfo.value.code == 1
        assert "uncommitted or untracked" in capsys.readouterr().out

    def test_no_commit_history_returns_none(self, tmp_path, monkeypatch):
        def fakeRun(cmd, **kwargs):
            if "rev-parse" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="true\n")
            if "status" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="")
            if "log" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="")
            raise AssertionError(f"unexpected command {cmd}")

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        assert misc.configGitInfo(str(tmp_path / "brb.ini")) is None

    def test_uses_custom_git_binary(self, tmp_path, monkeypatch):
        seen = []

        def fakeRun(cmd, **kwargs):
            seen.append(cmd[0])
            if "rev-parse" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="true\n")
            if "status" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="")
            if "log" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="abc\n")
            raise AssertionError(f"unexpected command {cmd}")

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        misc.configGitInfo(str(tmp_path / "brb.ini"), gitBin="/usr/local/bin/git")
        assert all(c == "/usr/local/bin/git" for c in seen)

    def test_status_timeout_warns_and_returns_none(self, tmp_path, monkeypatch, capsys):
        def fakeRun(cmd, **kwargs):
            if "rev-parse" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="true\n")
            if "status" in cmd:
                raise sp.TimeoutExpired(cmd, 5)
            raise AssertionError(f"unexpected command {cmd}")

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        assert misc.configGitInfo(str(tmp_path / "brb.ini")) is None
        err = capsys.readouterr().err
        assert "could not run git status" in err

    def test_status_calledprocesserror_warns_and_returns_none(
        self, tmp_path, monkeypatch, capsys
    ):
        def fakeRun(cmd, **kwargs):
            if "rev-parse" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="true\n")
            if "status" in cmd:
                raise sp.CalledProcessError(128, cmd)
            raise AssertionError(f"unexpected command {cmd}")

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        assert misc.configGitInfo(str(tmp_path / "brb.ini")) is None
        assert "could not run git status" in capsys.readouterr().err

    def test_git_binary_vanishes_after_probe_returns_none(self, tmp_path, monkeypatch):
        def fakeRun(cmd, **kwargs):
            if "rev-parse" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="true\n")
            raise FileNotFoundError

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        assert misc.configGitInfo(str(tmp_path / "brb.ini")) is None

    def test_log_timeout_warns_and_returns_none(self, tmp_path, monkeypatch, capsys):
        def fakeRun(cmd, **kwargs):
            if "rev-parse" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="true\n")
            if "status" in cmd:
                return sp.CompletedProcess(cmd, 0, stdout="")
            if "log" in cmd:
                raise sp.TimeoutExpired(cmd, 5)
            raise AssertionError(f"unexpected command {cmd}")

        monkeypatch.setattr(misc.sp, "run", fakeRun)
        assert misc.configGitInfo(str(tmp_path / "brb.ini")) is None
        assert "could not run git log" in capsys.readouterr().err

    def test_dirty_check_still_fails_loudly(self, tmp_path, monkeypatch, capsys):
        """A git failure must not silence a real dirty-config abort."""
        monkeypatch.setattr(
            misc.sp,
            "run",
            lambda cmd, **kw: sp.CompletedProcess(
                cmd,
                0,
                stdout="true\n" if "rev-parse" in cmd else " M brb.ini\n",
            ),
        )
        with pytest.raises(SystemExit) as excinfo:
            misc.configGitInfo(str(tmp_path / "brb.ini"))
        assert excinfo.value.code == 1
        assert "uncommitted or untracked" in capsys.readouterr().out


class TestPacifier:
    def test_removes_umlauts(self):
        assert misc.pacifier("Müller") == "Muller"

    def test_removes_spaces(self):
        assert misc.pacifier("John Doe") == "JohnDoe"

    def test_removes_apostrophes(self):
        assert misc.pacifier("O'Brien") == "OBrien"

    def test_plain_ascii_is_unchanged(self):
        assert misc.pacifier("PlainName") == "PlainName"


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
