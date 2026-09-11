import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from dominate.tags import br, div, html
from tabulate import tabulate

from BRB.logger import log
from BRB.misc import getVersion


def errorEmail(config, errTuple, msg):
    configCommit = config.get("Options", "configCommit", fallback="")
    msg = MIMEText(
        msg
        + f"\nError type: {errTuple[0]}\nError value: {errTuple[1]}\n{errTuple[2]}\n"
        + (f"\nConfig file commit: {configCommit}\n" if configCommit else "")
    )
    gitBin = config.get("software", "git", fallback="git")
    msg["Subject"] = f"[BigRedButton {getVersion('BRB', gitBin)}] Error"
    msg["From"] = config.get("Email", "fromAddress")
    msg["To"] = config.get("Email", "errorTo")

    s = smtplib.SMTP(config.get("Email", "host"))
    s.send_message(msg)
    s.quit()


def finishedEmail(config, msg):
    mailer = MIMEMultipart("alternative")
    gitBin = config.get("software", "git", fallback="git")
    mailer["Subject"] = (
        f"[BigRedButton {getVersion('BRB', gitBin)}] "
        f"{config.get('Options', 'runID')} processed"
    )
    mailer["From"] = config.get("Email", "fromAddress")

    # Create the table head
    _html = html()
    # Default recipient is finishedTo (bioinfocore)
    recipient = config.get("Email", "finishedTo")
    # Inform deepseq too if we have a sambaUpdate:
    if any(i[6] for i in msg):
        log.info("At least one sambaUpdate true in msg")
        # Only inform deepseq if no workflow failed, and no group was
        # skipped (a flowcell with a skipped group -- e.g. owned by a live
        # PID from another process -- was not actually fully processed this
        # pass, so deepSeq should not be told it's ready).
        statuses = [i[4] for i in msg]
        if statuses.count("FAILED") == 0 and not any(
            s.startswith("SKIPPED") for s in statuses
        ):
            recipient = config.get("Email", "deepSeq")
            _html.add(
                div(
                    f"Post-processing is ready, Samba drive is updated for {[i[6] for i in msg].count(True)} project(s).",
                    br(),
                )
            )

    mailer["To"] = recipient
    # Table
    tabHead = [
        "Project",
        "organism",
        "libraryType",
        "workflow",
        "workflow_status",
        "parkour_status",
        "sambaUpdate",
        "reruns",
    ]
    configCommit = config.get("Options", "configCommit", fallback="")
    message = (
        _html.render()
        + "\n\n"
        + tabulate(msg, tabHead, tablefmt="html", disable_numparse=True)
        + (f"\n\n<p>Config file commit: {configCommit}</p>" if configCommit else "")
    )

    email = MIMEText(message, "html")
    mailer.attach(email)

    s = smtplib.SMTP(config.get("Email", "host"))

    s.sendmail(
        config.get("Email", "fromAddress"), recipient.split(","), mailer.as_string()
    )
    s.quit()
