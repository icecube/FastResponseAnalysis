r'''Script to initiate fast reponse
    followups. Arguments can be passed
    at the beginning of execution or as 
    user I/O 

    Author: Alex Pizzuto
    Date: April 9, 2019
    '''

import argparse
import subprocess
import warnings
import fast_response.web_utils as web_utils
import logging as log
from fast_response.ExternalFollowup import ExternalFollowup, ExternalSkymapFollowup
warnings.filterwarnings("ignore")
log.basicConfig(level=log.ERROR)

parser = argparse.ArgumentParser(description='Fast Response Analysis')
parser.add_argument('--name', type=str, default=None,
                    help='Name of the source')
parser.add_argument('--skymap', type=str, default=None,
                    help='Source localization PDF path')
parser.add_argument('--ra', default=None, type=float,
                    help='Right ascension (in degrees)')
parser.add_argument('--dec', default=None, type=float,
                    help='Declination (in degrees)')
parser.add_argument('--start', type=str, required=True,
                    help="Start of the onsource time window relative "
                    "to the trigger in seconds")
parser.add_argument('--stop', type=str, required=True,
                    help="End of the onsource time window relative "
                    "to the trigger in seconds")
parser.add_argument('--extension', type=float, default=None,
                    help="Source extension in degrees")
parser.add_argument('--skip-events', default=None,
                    type= lambda z:[ tuple(int(y) for y in x.split(':')) for x in z.split(',')],
                    help="list of events to exclude from this analysis. "
                    "such as HESE events that contributed to the trigger."
                    "Example --skip-events  127853:67093193,128290:6888376")
parser.add_argument('--ntrials', default=100, type=int,
                    help="Number of background trials to perform")
parser.add_argument('--document', default=False, action='store_true')
parser.add_argument('--nodiag', default=False, action='store_true')
args = parser.parse_args()

assert ((args.skymap is None) != (args.ra is None and args.dec is None)), 'Must pass exactly one option between RA, Dec or skymap location'

subprocess.call(['clear'])
try:
    import pyfiglet
    message = '*'*80
    message += '\n' + str(pyfiglet.figlet_format("Fast Response Followup")) + '\n'
    message += '*'*80
    print(message)
except:
    message = '*'*80
    message += '\n' + 'Fast Response Followup' + '\n'
    message += '*'*80
    print(message)

if args.skymap is not None:
    f = ExternalSkymapFollowup(
        args.name, args.skymap, args.start, args.stop, skipped=args.skip_events,
        save=True, extension=None
    )
else:
    f = ExternalFollowup(
        args.name, args.ra, args.dec, args.start, args.stop, extension=args.extension,
        skipped=args.skip_events, save=True
    )

print(str(f))

f.unblind_TS()
f.plot_ontime()
f.ns_scan()
f.find_coincident_events()
f.calc_pvalue(ntrials = args.ntrials)
if not args.nodiag:
    f.make_dNdE()
    f.plot_tsd()
    f.upper_limit()
results = f.save_results()

print(results)
print(results.keys())

f.generate_report()
# if args.document:
#     subprocess.call(['cp','-r',results['analysispath'],
#         '/home/apizzuto/public_html/FastResponse/webpage/output/{}'.format(results['analysisid'])])
#     web_utils.updateFastResponseWeb(results)
