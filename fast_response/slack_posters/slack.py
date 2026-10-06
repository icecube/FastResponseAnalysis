r"""
Updated slack poster class to use new slack apps
Jessie Thwaites, March 2025
"""

import os
import pandas as pd
from datetime import datetime
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError
from dotenv import load_dotenv
from icecube import realtime_tools
load_dotenv()

class slackbot(object):
    """ Slackbot for sending messages
    
    Parameters
    ----------
    channel: str
        Channel name to post message
    """

    def __init__(self, channel_name):
        self.slack_api_key = os.getenv('SLACK_API_KEY')
        self.botname = 'fra-alert-bot'
        self.icon='https://upload.wikimedia.org/wikipedia/commons/5/5d/Penguu_%28Adelie_Penguin%29.png'

        self.channel_id = realtime_tools.slackapp_tools.get_channel_id(self.channel_name)
        if channel_id is not None:
            self.channel_id = channel_id
            self.channel_name = channel_name

    def get_shifter_id(self):
        analysis_path = os.environ.get('FAST_RESPONSE_SCRIPTS')
        on_shift = ''
        try: 
            shifters = pd.read_csv(os.path.join(analysis_path,'../slack_posters/fra_shifters.csv'), 
                                parse_dates=[0,1])
            for i in shifters.index:
                if shifters['start'][i]<datetime.utcnow()<shifters['stop'][i]:
                    on_shift+='<@{}> '.format(shifters['slack_id'][i]) 
        except Exception as e: 
            print(f"Error finding shifters: {e}")

        return on_shift

    def post_block_to_slack(self, message_dict, title='Test', header=False):
        """
        Format and post a formatted block message to slack
        
        Parameters
        ----------
        message_dict: dict
            Dictionary of items to post on slack (format {field_name: field_value})
            Note: The field_name will be bolded in the message
        title: str
            Title of the message to use. Only included if header=True
        header: bool
            Option to include a title and bot name in the message
        """ 
        if header:
            slack_block = [{
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": "*{0}* reports: {1}".format(self.botname,title)}
                    }]
        else:
            slack_block = []
            
        additional_fields = []
        for key, value in message_dict.items():
            additional_fields.append(
                {
                    "type": "mrkdwn",
                    "text": "*{0}*\n{1}".format(key, value)
                }
            )
        slack_block.append({"type" : "section", "fields": additional_fields})

        client = WebClient(token=self.slack_api_key)
        if self.channel_name is None:
            print('No channel specified, posting to test_messaging')
            self.channel_name = 'test_messaging'
            self.channel_id = realtime_tools.slackapp_tools.get_channel_id(self.channel_name)
        if self.channel_id is None:
            # should never get here, but if it does... exit
            print('Channel ID not found, exiting')
            return

        try:    
            # Call the chat.postMessage method using the WebClient
            result = client.chat_postMessage(
                channel=self.channel_id,
                text='Block message',
                blocks=slack_block,
                username=self.botname,
                icon_url=self.icon,
            )
            print('Message posted successfully to {}'.format(self.channel_name))

        except SlackApiError as e:
            print(f"Error posting message: {e}")

    def post_short_msg(self, message):
        """
        Post a short message to slack
        
        Parameters
        ----------
        message: str
            String to be posted (allows links or newlines if needed)
        """ 

        client = WebClient(token=self.slack_api_key)
        if self.channel_name is None:
            print('No channel specified, posting to test_messaging')
            self.channel_name = 'test_messaging'
            self.channel_id = realtime_tools.slackapp_tools.get_channel_id(self.channel_name)
        if self.channel_id is None:
            # should never get here, but if it does... exit
            print('Channel ID not found, exiting')
            return
        
        try:    
            # Call the chat.postMessage method using the WebClient
            result = client.chat_postMessage(
                channel=self.channel_id,
                text=message,
                username=self.botname,
                icon_url=self.icon,
            )
            print('Message posted successfully to {}'.format(self.channel_name))
            
        except SlackApiError as e:
            print(f"Error posting message: {e}")
 
    def post_file_to_slack(self, title = "results file", file_name=None):
        """
        Post a file to a slack channel

        Parameters
        ----------
        title: str
            title of the file, included in the post
        file_name: str
            path to the file being sent (required)
        """
        if file_name is None:
            print('No file specified, exiting')
            return
        
        client = WebClient(token=self.slack_api_key)
        if self.channel_name is None:
            print('No channel specified, posting to test_messaging')
            self.channel_name = 'test_messaging'
            self.channel_id = realtime_tools.slackapp_tools.get_channel_id(self.channel_name)
        if self.channel_id is None:
            # should never get here, but if it does... exit
            print('Channel ID not found, exiting')
            return
        
        comment = "*{0}* reports: {1}".format(self.botname,title)
        try:
            # Call the chat.postMessage method using the WebClient
            result = client.files_upload_v2(
                file = file_name,
                channel=self.channel_id,  
                title=file_name,
                initial_comment=comment,
                username=self.botname,
                icon_url = self.icon,
            )
            print('File posted successfully to {}'.format(self.channel_name))
            
        except SlackApiError as e:
            print(f"Error posting file: {e}")
