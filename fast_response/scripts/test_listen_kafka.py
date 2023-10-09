#!/usr/bin/env python

import logging
from gcn_kafka import Consumer
from icecube import realtime_tools
import json
import argparse

## none of these bools work yet. just listening to the real one
#parser = argparse.ArgumentParser(description='test listener for icecube kafka fra/llama results')
#parser.add_argument('--test_domain', type=bool, default=False,
#                        help='bool to use test.gcn.nasa.gov (default False)')
#parser.add_argument('--test_topic', type=bool, default=True,
#                        help='listen to gcn.notices.icecube.TEST.lvk_nu_track_search')
#args = parser.parse_args()


with open('/home/jthwaites/private/tokens/kafka_token.txt') as f:
    client_id = f.readline().rstrip('\n')
    client_secret = f.readline().rstrip('\n')

#domain = 'test.gcn.nasa.gov'
domain = 'gcn.nasa.gov'

consumer = Consumer(client_id=client_id,
                    client_secret=client_secret,
                    domain=domain)

#topic = 'gcn.notices.icecube.test.lvk_nu_track_search'
topic = 'gcn.notices.icecube.lvk_nu_track_search'

consumer.subscribe([topic])

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.warning("checking for {}, connecting to GCN".format(topic))

while True:
    for message in consumer.consume(timeout=1):
        if message.error():
            print(message.error())
            continue
        value = message.value()
        
        alert_dict = json.loads(value.decode('utf-8'))
        print(json.dumps(alert_dict, indent=2))
    
