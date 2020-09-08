import sys
sys.path.append('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/time_integrated_scripts/')
from steady_sensitivity_fits import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

output_path = '/data/user/apizzuto/fast_response_skylab/alert_event_followup/alert_skymaps/'

skymap_files = glob('/data/ana/realtime/alert_catalog_v2/fits_files/Run*.fits.gz')
l_ind = skymap_files[0].find("Run")
r_ind = skymap_files[0].find("_nside")

for ind in range(188, len(skymap_files)):
    title = skymap_files[ind][l_ind:r_ind].replace('_', '_Event_').replace('Run', 'Run_')
    plot_zoom(ind, LLH=True)
    zoom_str = 'zoom_LLH'
    plt.savefig(output_path + '{}_{}.png'.format(title, zoom_str), bbox_inches='tight')
    plt.close()

    #plot_zoom(ind, LLH=False)
    #zoom_str = 'probs'
    #plt.savefig(output_path + '{}_{}.png'.format(title, zoom_str), bbox_inches='tight')
    #plt.close()

    plot_skymap(ind)
    zoom_str = 'allsky'
    plt.savefig(output_path + '{}_{}.png'.format(title, zoom_str), bbox_inches='tight')
    plt.close()    
