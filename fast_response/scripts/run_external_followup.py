r'''Script to initiate fast reponse
    followups.

    Author: Alex Pizzuto
    Date: 2021
    '''

from fast_response.ExternalFollowup import ExternalFollowup, ExternalSkymapFollowup
import argparse
import subprocess
import warnings
import fast_response.web_utils as web_utils
import logging as log
import pyfiglet

def run_analysis(args):
    assert ((args.skymap is None) != (args.ra is None and args.dec is None)), 'Must pass exactly one option between RA, Dec or skymap location'

    #subprocess.call(['clear'])
    message = '*'*80
    message += '\n' + str(pyfiglet.figlet_format("Fast Response Followup")) + '\n'
    message += '*'*80
    print(message)

    if args.skymap is not None:
        f = ExternalSkymapFollowup(
            args.name, args.skymap, args.start,
            args.stop, skipped=args.skip_events,
            save=True, extension=None
        )
    else:
        f = ExternalFollowup(
            args.name, args.ra, args.dec, args.start,
            args.stop, extension=args.extension,
            skipped=args.skip_events, save=True
        )

    print(str(f))

    f.unblind_TS()
    f.plot_ontime()
    f.ns_scan()
    f.find_coincident_events()
    f.calc_pvalue(ntrials = args.ntrials)
    f.make_dNdE()
    f.plot_tsd()
    f.upper_limit(n_per_sig = args.n_per_sig)
    results = f.save_results()

    print("Beginning to make the report")

    f.generate_report()
    # if args.document:
    #     subprocess.call(['cp','-r',results['analysispath'],
    #         '/home/apizzuto/public_html/FastResponse/webpage/output/{}'.format(results['analysisid'])])
    #     web_utils.updateFastResponseWeb(results)

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    log.basicConfig(level=log.ERROR)

    parser = argparse.ArgumentParser(description='Fast Response Analysis')
    parser.add_argument('--name', type=str, default="Fast Response",
                        help='Name of the source')
    parser.add_argument('--skymap', type=str, default=None,
                        help='Source localization PDF path')
    parser.add_argument('--ra', default=None, type=float,
                        help='Right ascension (in degrees)')
    parser.add_argument('--dec', default=None, type=float,
                        help='Declination (in degrees)')
    parser.add_argument('--start', type=str, required=True,
                        help="Start time of the analysis in ISO format")
    parser.add_argument('--stop', type=str, required=True,
                        help="End time of the analysis in ISO format")
    parser.add_argument('--extension', type=float, default=None,
                        help="Source extension in degrees")
    parser.add_argument('--skip-events', default=None,
                        type= lambda z:[ tuple(int(y) for y in x.split(':')) for x in z.split(',')],
                        help="Event to exclude from the analyses, eg."
                        "Example --skip-events=127853:67093193")
    parser.add_argument('--ntrials', default=100, type=int,
                        help="Number of background trials to perform")
    parser.add_argument('--n_per_sig', default=100, type=int,
                        help="Number of signal trials per injected ns to perform")
    parser.add_argument('--document', default=False, action='store_true')
    args = parser.parse_args()

    run_analysis(args)
