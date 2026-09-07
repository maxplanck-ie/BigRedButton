import configparser
import threading
import time

import BRB.ET


def parkourConfig(tmp_path):
    config = configparser.ConfigParser()
    config["Paths"] = {"baseData": str(tmp_path)}
    config["Options"] = {"runID": "210608_A00931_0309_BHCCMWDRXY"}
    config["Parkour"] = {
        "ResultsURL": "https://parkour.example.org/api/results/",
        "user": "u",
        "password": "p",
        "cert": "",
    }
    return config


class TestSendToParkourSerialisation:
    def test_concurrent_posts_are_serialised(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        state = {"inFlight": 0, "maxInFlight": 0, "posts": 0}
        guard = threading.Lock()

        def fakePost(url, auth=None, data=None, verify=None):
            with guard:
                state["inFlight"] += 1
                state["posts"] += 1
                state["maxInFlight"] = max(state["maxInFlight"], state["inFlight"])
            time.sleep(0.05)
            with guard:
                state["inFlight"] -= 1
            return "RESPONSE"

        monkeypatch.setattr(BRB.ET.requests, "post", fakePost)

        threads = [
            threading.Thread(
                target=BRB.ET.sendToParkour, args=(config, [{"barcode": str(i)}])
            )
            for i in range(4)
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=10)

        assert state["posts"] == 4
        assert state["maxInFlight"] == 1

    def test_return_value_unchanged(self, tmp_path, monkeypatch):
        config = parkourConfig(tmp_path)
        monkeypatch.setattr(BRB.ET.requests, "post", lambda *a, **k: "RESPONSE")
        assert BRB.ET.sendToParkour(config, [{"barcode": "1"}]) == "RESPONSE"


class TestGetOffSpeciesRate:
    def _writeReport(self, d, lines):
        (d / "sample.rep").write_text("\n".join(lines) + "\n")

    def test_known_label_matches_its_kraken_group(self, tmp_path):
        self._writeReport(
            tmp_path,
            ["97.00\t123\t45\tG1\t0\tflygrp", "1.00\t1\t1\tG1\t0\tother"],
        )
        rate = BRB.ET.getOffSpeciesRate(str(tmp_path), "drosophila")
        assert rate == 1 - 97.00 / 100

    def test_unknown_label_returns_zero_and_does_not_raise(self, tmp_path):
        self._writeReport(tmp_path, ["97.00\t123\t45\tG1\t0\tflygrp"])
        assert BRB.ET.getOffSpeciesRate(str(tmp_path), "ricefish") == 0
        assert BRB.ET.getOffSpeciesRate(str(tmp_path), None) == 0

    def test_known_label_whose_group_is_absent_from_report_returns_zero(self, tmp_path):
        # regression: org_map[org_label] never matches a line in the report,
        # `off` used to stay unbound and raise UnboundLocalError instead of
        # returning a sane default.
        self._writeReport(tmp_path, ["100.00\t123\t45\tG1\t0\tsomeothergrp"])
        assert BRB.ET.getOffSpeciesRate(str(tmp_path), "drosophila") == 0
