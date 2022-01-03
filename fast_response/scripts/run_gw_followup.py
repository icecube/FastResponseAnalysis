import argparse, subprocess
from astropy.time import Time
import pyfiglet

from fast_response.GWFollowup import GWFollowup

def run_gw_followup(name, time, skymap):
    subprocess.call(['clear'])
    message = '*'*80
    message += '\n' + str(pyfiglet.figlet_format("GW Followup")) + '\n'
    message += '*'*80
    print(message)

    gw_time = Time(time, format='mjd')

    delta_t = 1000.
    start_time = gw_time - (delta_t / 86400. / 2.)
    stop_time = gw_time + (delta_t / 86400. / 2.)
    start = start_time.iso
    stop = stop_time.iso

    name = name.replace('_', ' ')
    name = name + 'new code framework'

    f = GWFollowup(name, skymap, start, stop)
    f.unblind_TS()
    f.plot_ontime()
    f.calc_pvalue()
    f.make_dNdE()
    f.plot_tsd()
    f.upper_limit()
    f.find_coincident_events()
    f.per_event_pvalue()
    results = f.save_results()
    f.generate_report()

    f.write_circular()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='GW Followup')
    parser.add_argument('--skymap', type=str, default=None,
                        help='path to skymap')
    parser.add_argument('--time', type=float, default=None,
                        help='Time of the GW (mjd)')
    parser.add_argument('--name', type=str,
                        default="GW Followup")
    args = parser.parse_args()

    run_gw_followup(args.name, args.time, args.skymap)