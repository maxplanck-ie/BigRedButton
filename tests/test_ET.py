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
