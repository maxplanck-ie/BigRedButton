import smtplib
from email.mime.text import MIMEText
from importlib.metadata import version

def errorEmail(config, errTuple, msg) :
    msg = MIMEText(msg + "\nError type: %s\nError value: %s\n%s\n" % (errTuple[0], errTuple[1], errTuple[2]))
    msg['Subject'] = f"[BigRedButton {version("BRB")}] Error"
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","errorTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()


def finishedEmail(config, msg) :
    message = "Flow cell: {}\n".format(config.get("Options","runID"))
    message += msg

    msg = MIMEText(message)
    msg['Subject'] = f"[BigRedButton {version("BRB")}] {config.get("Options","runID")} processed"
    msg['From'] = config.get("Email","fromAddress")
    msg['To'] = config.get("Email","finishedTo")

    s = smtplib.SMTP(config.get("Email","host"))
    s.send_message(msg)
    s.quit()
