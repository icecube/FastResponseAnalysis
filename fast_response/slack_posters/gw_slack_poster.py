'''
Script to make a slackbot and post 
Slack notifications for GW alerts
11/18 - currently for mocks
Nov 2022, Jessie Thwaites
'''

with open('mock_gw_slackbot.txt') as f:
    channel = f.readline()
    webhook = f.readline()
    bot_name = f.readline()

import gcn
@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE)

def process_gcn(payload, root):
    
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    # Read trigger time of event
    eventtime = root.find('.//ISOTime').text
    AlertTime = Time(eventtime, format='isot').iso

    name = root.attrib['ivorn'].split('#')[1]
    if 'MS22' in name: mock=True
    else: mock=False

    bot = slackbot(channel, bot_name, webhook)
    if mock: 
        bot.send_message(f'Mock GW alert found: {name}, {AlertTime} UTC', 'chaos')
    else: 
        bot.send_message(f'GW alert found: {name}, {AlertTime} UTC', 'chaos')
    print('\n \n')

if __name__ == '__main__':
    import sys
    from slack import slackbot
    from astropy.time import Time
    import argparse

    #parser = argparse.ArgumentParser(description='FRA GW followup slackbot')
    #parser.add_argument('--testing', action='store_true', default=False,
    #                    help='run a test message')
    #args = parser.parse_args()
    print('Listening for GCNs . . .')
    gcn.listen(handler=process_gcn)

