r"""

@file :   slack.py
@author:  Josh Wood
@date:    3 Oct 2014
@brief:   class for generating a slackbot and having it send messages to slack.

"""

import os

class slackbot(object):

    def __init__(self, channel, name, url):
        r""" Constructor

        Parameters
        ----------
        channel : str
            Channel name to send to
        name : str
            Name of author who is sending messages
            (Note: doesn't need to be a user account)
        url : str
            Url path for webhook.

        NOTE: To create URL for webhook, go to your Slack desktop app.
        Under the 'IceCube' drop down menu choose 'Customize Slack'.
        This brings you to a webpage with a menu in the upper left.
        Click Menu --> Configure Apps --> Custom Integrations
        --> Incoming WebHooks. Click the green "add configuration"
        button on the left of the Incoming WebHooks page. Search
        for a channel name to post to and click Add.
        """
        self.name = name
        self.channel = channel
        self.url = url

    def send_message(self, message, emoji):
        r""" Sends message to slack channel specified
             by self.channel under the name self.name.
             Also posts an emoji next to first message
             in the thread.

        Parameters
        ----------
        message : str
            One line message. alerts should never be longer.
        emoji : str
            Emoji icon to appear next to message
        """

        # command needed to post data to the web
        cmd = "curl -X POST --data-urlencode"

        # payload of the slack message
        pay = (('payload={"channel": "%s", "username": "%s", ' +
                '"text": "%s", "icon_emoji": ":%s:"}') %
               (self.channel, self.name, message, emoji))

        # send message with full system command
        syscmd = cmd + " \'" + pay + "\' " + self.url
        os.system(syscmd)

