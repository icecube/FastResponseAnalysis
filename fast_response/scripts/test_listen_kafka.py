#!/usr/bin/env python

import logging
from gcn_kafka import Consumer
from icecube import realtime_tools
import json, argparse

parser = argparse.ArgumentParser(description='test listener for icecube kafka fra/llama results')
parser.add_argument('--test_domain', action='store_true', default=False,
                    help='bool to use test IceCube stream (default False)')
parser.add_argument('--use_prod', action='store_true', default=False,
                    help='use production token rather than read-only')
parser.add_argument('--save_out', action='store_true', default=False,
                    help='save the packet as a test file')
parser.add_argument('--classic', action='store_true', default=False,
                    help='listen to the voevent streams from gcn rather than the kafka')
args = parser.parse_args()

logger = logging.getLogger()
logger.setLevel(logging.INFO)

if args.use_prod:
    token = '/home/jthwaites/private/tokens/real_icecube_kafka_prod.txt'
else:
    token = '/home/jthwaites/private/tokens/kafka_token.txt'

with open(token) as f:
    client_id = f.readline().rstrip('\n')
    client_secret = f.readline().rstrip('\n')

if args.test_domain:
    domain = 'test.gcn.nasa.gov'
else:
    domain = 'gcn.nasa.gov'

config = {'broker.address.family': 'v4', 
          'log_level': 0}
consumer = Consumer(client_id=client_id,
                    client_secret=client_secret,
                    domain=domain,
                    config = config
                   )

# choose topics to listen to
if args.classic:
    #no lvk_nu_track_search classic version
    topics = ['gcn.classic.voevent.ICECUBE_ASTROTRACK_BRONZE',
              'gcn.classic.voevent.ICECUBE_ASTROTRACK_GOLD',
              'gcn.classic.voevent.ICECUBE_CASCADE']
else:
    topics =['gcn.notices.icecube.gold_bronze_track_alerts',
             'gcn.notices.icecube.test.gold_bronze_track_alerts',
             'gcn.notices.icecube.lvk_nu_track_search',
             'igwn.gwalert']

consumer.subscribe(topics)
logger.warning("checking for alerts, connecting to GCN")#.format(topic))

while True:
    for message in consumer.consume(timeout=1):
        if message.error():
            print(message.error())
            continue
        value = message.value()
        logger.warning('Found GCN on topic {}'.format(message.topic()))

        try:
            if args.classic:
                alert_xml = value.decode('utf-8') #.encode('ascii')
                print(alert_xml)
                if args.save_out:
                    with open('test_alert.xml', 'w') as f:
                        f.write(alert_xml)
            else:
                alert_dict = json.loads(value.decode('utf-8'))
                if args.save_out:
                    with open('test_alert.json', 'w') as f:
                        json.dump(alert_dict, f)
                if 'skymap' in alert_dict:
                    alert_dict.pop('skymap', None)
                print(json.dumps(alert_dict, indent=2))

        except Exception as e:
            print('Caught exception when decoding: ', e, '\nSkipping...')
            continue
        
