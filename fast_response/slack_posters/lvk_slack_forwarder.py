#!/usr/bin/env python

import logging
from datetime import datetime
import socket
import gcn
import requests
import healpy as hp
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import io, time
import urllib.request, urllib.error, urllib.parse
import argparse

from icecube import realtime_tools

# Function to call every time a GCN is received.
# Run only for notices of type
# LVC_PRELIMINARY, LVC_INITIAL, or LVC_UPDATE.
@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_EARLY_WARNING,
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE,
    gcn.notice_types.LVC_RETRACTION)

def process_gcn(payload, root):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.info("lvc_slack_forwarder processing GCN")

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}
    
    # Where to post:
    heartbeat=False
    if 'Significant' in params.keys():
        if int(params['Significant'])==0:
            #subthreshold alerts
            heartbeat=True

    if realtime_tools.config.TESTING:
        slack_channel = ['#test_messaging']
    elif root.attrib['role'] == 'observation' and not heartbeat:
        slack_channel = ['#alerts', '#alerts-heartbeat', '#gwnu-heartbeat','#fra-shifting']
    elif root.attrib['role'] == 'observation':
        slack_channel = ['#alerts-heartbeat', '#gwnu-heartbeat']
    else:
        slack_channel = ['#gw-mock-heartbeat']

    #If retracted - most parameters will be missing. Post link to GraceDB and event name only
    if params['AlertType'] =='Retraction':
        slack_message = {
            "icon_emoji": ":gw:",
            "text" : "Event {0} Retracted. See URL: {1}.".format(params['GraceID'], params['EventPage']),
        } 
        for channel in slack_channel:
            realtime_tools.messaging.to_slack(channel=channel,
                                      username="lvc-gcn-bot",
                                      content=slack_message["text"],
                                      **slack_message)
        return
    
    # Respond only to 'CBC' events. Change 'CBC' to "Burst'
    # to respond to only unmodeled burst events.
    #if params['Group'] != 'CBC':
    #    return

    lvc_params = {}
    lvc_params["Alert Time"] = root.find("./WhereWhen//ISOTime").text
    lvc_params["GraceDB ID"] = params['GraceID']
    lvc_params["Alert Type"] = params['AlertType']
    lvc_params["GraceDB URL"] = params['EventPage']
    lvc_params["False alarm rate"] = '{:.4g}'.format(float(params['FAR']))
    lvc_params["Pipeline"] = params['Pipeline']
    ## Some of these params are pipeline specfic, check they're there.
    if 'BBH' in params:
        lvc_params["BBH Prob"] = round(float(params['BBH']),4)
    if 'BNS' in params:    
        lvc_params["BNS Prob"] = round(float(params['BNS']),4)
    if 'NSBH' in params:
        lvc_params["NSBH Prob"] = round(float(params['NSBH']),4)
    if 'Noise' in params:
        lvc_params['Terrestrial Prob'] = round(float(params['Noise']),4)
    if 'Terrestrial' in params:
        lvc_params["Terrestrial Prob"] = round(float(params['Terrestrial']),4)
    if 'Significant' in params:
        lvc_params['Significant'] = 'True' if int(params['Significant'])==1 else 'False'
    
    #properties to tell if there is a neutron star
    if 'HasRemnant' in params:
        lvc_params["HasRemnant Prob"] = round(float(params['HasRemnant']),4)
    if 'HasNS' in params:
        lvc_params['HasNS Prob'] = round(float(params['HasNS']),4)
    if 'HasMassGap' in params:
        lvc_params['MassGap Prob'] = round(float(params['HasMassGap']),4)
    
    if 'CentralFreq' in params:
        lvc_params["CWB Central Freq"] = round(float(params['CentralFreq']),4)
    if 'Duration' in params:
        lvc_params["CWB Duration"] = round(float(params['Duration']),4)    
        
    if 'skymap_fits' in params:
        # Read the HEALPix sky map and the FITS header.
        skymap_link = params['skymap_fits']
        #if we don't wait, the old format isn't uploaded
        time.sleep(6.)
        if 'multiorder' in skymap_link:
            try:
                bayestar_map=skymap_link.replace('multiorder.','').split(',')
                if len(bayestar_map)==1:
                    new_map=bayestar_map[0]+'.gz'
                else: 
                    new_map=bayestar_map[0]+'.gz,'+ bayestar_map[1]
                ret = requests.head(new_map)
                assert ret.status_code == 200
                skymap_link=new_map
            except:
                logger.warning('Failed to download skymap in correct format')

        try:
            skymap, header = hp.read_map(skymap_link, h=True, verbose=False)
            header = dict(header)
            if 'DISTMEAN' in header:
                lvc_params["Src Distance"] = str(round(float(header['DISTMEAN']),3)) + \
                                             '+/-' + str(round(float(header['DISTSTD']),3))
                if 'INSTRUME' in header:
                    lvc_params["Instruments"] = header['INSTRUME']
        except Exception as e: #urllib.error.HTTPError:
            #logger.error("skymap_fits file not found in GraceDB!")
            logger.error("error opening skymap FITS file")
            skymap = None
            header = None

    slack_fields_start = []
    slack_fields_end = []  # Attach longer info to end of message
    for key,value in list(lvc_params.items()):
        if key in ['GraceDB URL']:
            entry = {
                "title": key,
                "value": value,
                "short": False
            }
            slack_fields_end.append(entry)
        else:
            entry = {
                "title": key,
                "value": value,
                "short": True
            }
            slack_fields_start.append(entry)
    slack_fields = slack_fields_start + slack_fields_end
    obs_type = params['AlertType']

    if root.attrib['role'] != 'observation':
        text = "Mock LVK alert received: Type: {0}".format(obs_type)
    else: 
        text = "LVK alert received: Type: {0}".format(obs_type)
    slack_message = {
        "icon_emoji": ":gw:",
        "text" : text,
        "color":"danger"
    }
    epoch = datetime.utcfromtimestamp(0)
    ts = int((datetime.utcnow() - epoch).total_seconds())
    slack_attach = [
        {
            "title": "Alert Information",
            "color": "danger",
            "mrkdwn_in": ["fields"],
            "fields": slack_fields,
            "footer": socket.gethostname(),
            "ts": ts,
        },
    ]
    slack_message["attachments"] = slack_attach

    # # if not Initial or Preliminary, just a brief post:
    if obs_type !='Preliminary' and obs_type != 'Initial' and obs_type != 'EarlyWarning':
        slack_message = {
            "icon_emoji": ":gw:",
            "text" : "LVK alert info received: Type: {0}. See URL: {1}.".format(obs_type,
                                                                                lvc_params["GraceDB URL"]),
        } 

    for channel in slack_channel:
        realtime_tools.messaging.to_slack(channel=channel,
                                      username="lvc-gcn-bot",
                                      content=slack_message["text"],
                                      **slack_message)
    
    '''#if there is an NS
    has_ns=False
    if ('BNS Prob' and 'NSBH Prob') in lvc_params.keys():
        if lvc_params["BNS Prob"]+lvc_params["NSBH Prob"]>0.5:
            has_ns=True
    elif ('HasRemnant Prob' or 'HasNS Prob') in lvc_params.keys():
        if lvc_params['HasRemnant Prob'] > 0.5:
            has_ns=True
        elif lvc_params['HasNS Prob'] > 0.5:
            has_ns=True
    
    if has_ns == True and not heartbeat:
        if obs_type == 'Initial' or obs_type =='Preliminary':
            realtime_tools.messaging.to_slack(channel='#alerts',
                                              username="lvc-gcn-bot",
                                              content=slack_message["text"],
                                              **slack_message)
            #print('Prob of having NS > 50%')
    '''

    # Make skymap, save to memfile
    # only post if Preliminary (or initial)
    if obs_type == 'Preliminary' or obs_type == 'Initial' or obs_type=='EarlyWarning':
        if skymap is not None:
            cmap = plt.cm.YlOrBr
            cmap.set_under('w')
            hp.mollview(skymap, title="", cbar=True, notext=True, rot=180, hold=False, cmap=cmap)
            hp.graticule()
            # add some labels
            plt.text(2.0,0., r"$0^\circ$", ha="left", va="center")
            plt.text(1.9,0.45, r"$30^\circ$", ha="left", va="center")
            plt.text(1.4,0.8, r"$60^\circ$", ha="left", va="center")
            plt.text(1.9,-0.45, r"$-30^\circ$", ha="left", va="center")
            plt.text(1.4,-0.8, r"$-60^\circ$", ha="left", va="center")
            plt.text(2.0, -0.15, r"$0\,\mathrm{h}$", ha="center", va="center")
            plt.text(1.333, -0.15, r"$4\,\mathrm{h}$", ha="center", va="center")
            plt.text(.666, -0.15, r"$8\,\mathrm{h}$", ha="center", va="center")
            plt.text(0.0, -0.15, r"$12\,\mathrm{h}$", ha="center", va="center")
            plt.text(-.666, -0.15, r"$16\,\mathrm{h}$", ha="center", va="center")
            plt.text(-1.333, -0.15, r"$20\,\mathrm{h}$", ha="center", va="center")
            plt.text(-2.0, -0.15, r"$0\,\mathrm{h}$", ha="center", va="center")
            #plt.savefig('./test_map.png', format='png', dpi=150)
            memfile = io.BytesIO()
            plt.savefig(memfile, format = 'png', dpi = 150)
        
            with open('/cvmfs/icecube.opensciencegrid.org/users/jthwaites/tokens/gw_token.txt') as f:
                my_key = f.readline()

            #post to slack
            for channel in slack_channel:
                memfile.seek(0)
                response = requests.post('https://slack.com/api/files.upload',
                                     timeout=60,
                                     params={'token': my_key},
                                     data={'filename':'lvk_skymap.png',
                                           'title': 'LVK GraceDB Skymap',
                                           'channels': channel},
                                     files={'file': memfile}
                                     )
                if response.ok is True:
                    logger.info("LVK skymap posted OK to {}".format(channel))
                else:
                    logger.error("Error posting skymap to {}!".format(channel))
                
        
# create the daemon functionality, this one listens to EVERYTHING
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='FRA GW followup')
    parser.add_argument('--testing', action='store_true', default=False,
                        help='Run on live GCNs')
    args = parser.parse_args()

    if not args.testing:
        # Listen for GCNs until the program is interrupted
        # (killed or interrupted with control-C).
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        logger.warning("lvc_slack_forwarder starting, connecting to GCN")
        
        gcn.listen(handler=process_gcn)

    else:
        # testing
        import lxml.etree
        import os, time
        
        #example files for all types live in /realtime_scripts/resources/test/
        files = ['S230524x-1-EarlyWarning.xml,0']
                 #'MS181101ab-1-Preliminary.xml',
                 #'S200308e-3-Retraction.xml']
   
        # excerise prelim, initial, update alerts:
        for file in files:
            src_path = os.environ['I3_SRC']
            src_path = '/data/user/jthwaites/realtime_updated'
            test_file = src_path+'/realtime_scripts/resources/test/'+file

            payload = open(test_file, 'rb').read()
            root = lxml.etree.fromstring(payload)
            root.attrib['role'] = 'test'

            process_gcn(payload, root)


