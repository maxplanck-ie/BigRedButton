import configparser
from unittest.mock import MagicMock, patch

from BRB import email


def make_config():
    config = configparser.ConfigParser()
    config["Options"] = {"runID": "run1"}
    config["Email"] = {
        "fromAddress": "brb@example.org",
        "finishedTo": "bioinfocore@example.org",
        "deepSeq": "deepseq@example.org",
        "errorTo": "errors@example.org",
        "host": "localhost",
    }
    return config


class TestErrorEmail:
    @patch("BRB.email.smtplib.SMTP")
    @patch("BRB.email.getVersion", return_value="v1.0.0")
    def test_sends_to_errorto_with_exception_details(self, mockGetVersion, mockSMTP):
        config = make_config()
        mockInstance = MagicMock()
        mockSMTP.return_value = mockInstance
        errTuple = (ValueError, ValueError("boom"), "traceback text")

        email.errorEmail(config, errTuple, "flowcell 210608_A00931 failed")

        assert mockInstance.send_message.called
        sentMsg = mockInstance.send_message.call_args[0][0]
        assert sentMsg["To"] == "errors@example.org"
        assert sentMsg["From"] == "brb@example.org"
        assert "v1.0.0" in sentMsg["Subject"]
        body = sentMsg.get_payload()
        assert "flowcell 210608_A00931 failed" in body
        assert "ValueError" in body
        assert "traceback text" in body
        assert mockInstance.quit.called

    @patch("BRB.email.smtplib.SMTP")
    @patch("BRB.email.getVersion", return_value="v1.0.0")
    def test_includes_config_commit_when_present(self, mockGetVersion, mockSMTP):
        config = make_config()
        config["Options"]["configCommit"] = "abc1234"
        mockInstance = MagicMock()
        mockSMTP.return_value = mockInstance
        errTuple = (RuntimeError, RuntimeError("x"), "tb")

        email.errorEmail(config, errTuple, "failure")

        sentMsg = mockInstance.send_message.call_args[0][0]
        assert "abc1234" in sentMsg.get_payload()

    @patch("BRB.email.smtplib.SMTP")
    @patch("BRB.email.getVersion", return_value="v1.0.0")
    def test_omits_config_commit_line_when_absent(self, mockGetVersion, mockSMTP):
        config = make_config()
        mockInstance = MagicMock()
        mockSMTP.return_value = mockInstance
        errTuple = (RuntimeError, RuntimeError("x"), "tb")

        email.errorEmail(config, errTuple, "failure")

        sentMsg = mockInstance.send_message.call_args[0][0]
        assert "Config file commit" not in sentMsg.get_payload()


class TestFinishedEmailSkippedSuppression:
    """
    A flowcell containing a SKIPPED-status entry was not actually fully
    processed this pass (see PushButton.runOneGroup's live-PID marker skip),
    so email.finishedEmail must not tell deepSeq it's ready to look at the
    Samba drive, even when some other group did update Samba.
    """

    @patch("BRB.email.smtplib.SMTP")
    def test_skipped_entry_suppresses_deepseq_notification(self, mockSMTP):
        config = make_config()
        mockInstance = MagicMock()
        mockSMTP.return_value = mockInstance
        msg = [
            [
                "1_A_Foo",
                "human",
                "ChIP-Seq",
                "DNA",
                "SKIPPED (owned by live PID)",
                "not updated",
                False,
                0,
            ],
            [
                "2_B_Bar",
                "mouse",
                "stranded mRNA-Seq",
                "RNA",
                "success",
                "PARKOUR_OK",
                True,
                0,
            ],
        ]

        email.finishedEmail(config, msg)

        assert mockInstance.sendmail.called
        recipients = mockInstance.sendmail.call_args[0][1]
        assert recipients == ["bioinfocore@example.org"]

    @patch("BRB.email.smtplib.SMTP")
    def test_no_skipped_or_failed_entries_still_notifies_deepseq(self, mockSMTP):
        config = make_config()
        mockInstance = MagicMock()
        mockSMTP.return_value = mockInstance
        msg = [
            [
                "2_B_Bar",
                "mouse",
                "stranded mRNA-Seq",
                "RNA",
                "success",
                "PARKOUR_OK",
                True,
                0,
            ],
        ]

        email.finishedEmail(config, msg)

        assert mockInstance.sendmail.called
        recipients = mockInstance.sendmail.call_args[0][1]
        assert recipients == ["deepseq@example.org"]
