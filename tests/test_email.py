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
        "host": "localhost",
    }
    return config


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
