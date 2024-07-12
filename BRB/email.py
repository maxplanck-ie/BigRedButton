import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from importlib.metadata import version
from dominate.tags import html, div, br
from tabulate import tabulate
from BRB.logger import log

def errorEmail(config, errTuple, msg) :
    msg = MIMEText(msg + "\nError type: %s\nError value: %s\n%s\n" % (errTuple[0], errTuple[1], errTuple[2]))
    msg['Subject'] = f"[BigRedButton {version("BRB")}] Error"
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","errorTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()


def finishedEmail(config, msg):
    mailer = MIMEMultipart('alternative')
    mailer['Subject'] = f"[BigRedButton {version("BRB")}] {config.get("Options","runID")} processed"
    mailer['From'] = config.get("Email","fromAddress")
    
    # Create the table head
    _html = html()
    # Default recipient is finishedTo (bioinfocore)
    recipient = config.get("Email","finishedTo")
    # Inform deepseq too if we have a sambaUpdate:
    if any([i[6] for i in msg]):
        log.info("At least one sambaUpdate true in msg")
        # Only inform deepseq if all workflows are succesfull (combat spam)
        if [i[4] for i in msg].count('success') == len(msg):
            recipient = config.get("Email","deepSeq")
            _html.add(div(
                "Post-processing is ready, deepSeq's sambda drive is updated for at least one project.",
                br()
            ))
    
    mailer['To'] = recipient
    # Table
    tabHead = ['Project', 'organism', 'libraryType', 'workflow', 'workflow_status', 'parkour_status', 'sambaUpdate']
    message = tabulate(
        msg, tabHead, tablefmt="html", disable_numparse=True
    )

    email = MIMEText(message, 'html')
    mailer.attach(email)
    
    s = smtplib.SMTP(config.get("Email", "host"))
    
    s.sendmail(
        config.get("Email","fromAddress"),
        recipient.split(','),
        mailer.as_string()
    )
    s.quit()