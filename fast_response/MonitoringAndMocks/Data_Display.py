#!/usr/bin/env python
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams['axes.titlesize'] = 18
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16
import matplotlib.pyplot as plt
import numpy as np
import pickle, glob, os
from datetime import date, timezone, datetime
from dateutil.relativedelta import *
# from statistics import median
from matplotlib.dates import DateFormatter
from astropy.time import Time
import pandas as pd
import subprocess, warnings
warnings.filterwarnings('ignore', module='astropy._erfa')

def dial_up(who="jessie"):
        cell_tower = "/home/jthwaites/private/"
        # subprocess.call([cell_tower+"make_call.py", f"--{who}=True", '--troubleshoot=True'])
        print('Calld')

path = os.environ.get('FAST_RESPONSE_OUTPUT')
out_put = '/data/user/jthwaites/o4-mocks/'

#Creating readable (and callable) files from ALL pickle files previously created in gw_gcn_listener
mock_files = sorted(glob.glob(path+'/PickledMocks/*MS*.pickle'), reverse=True)[:1000]

def sort_mocks(mock_files):
    event_dict = pd.DataFrame({"Trigger_Time": [], "GCN_Alert": [], "End_Time": [],
                "LVK_Latency": [], "IceCube_Latency": [], "Total_Latency": [], 
                "Name": [], "time_stamp": [], "first": [], 
                "last": [], "we_had_to_wait": []})
    ###these keys are saved in pickle files and just need to be loaded into dict
    select_keys = ("Trigger_Time", "GCN_Alert", "End_Time", "LVK_Latency",
                    "IceCube_Latency", "Total_Latency")
    ###create empty list and load each value into the list which is then loaded into the dict.
    loading_dict = {"Trigger_Time": [], "GCN_Alert": [], "End_Time": [],
                "LVK_Latency": [], "IceCube_Latency": [], "Total_Latency": [], 
                "Name": [], "time_stamp": [], "first": [], 
                "last": [], "we_had_to_wait": []}
    print('Beginning to load {} mocks'.format(len(mock_files)))
    for i in mock_files:
        mocks = pd.read_pickle(f'{i}')
        ###load data from pickle files to dataframe
        for key in select_keys:
            mock_key = "Ligo_Latency" if key == "LVK_Latency" else key
            loading_dict[key] = mocks[mock_key]

        ###separate name of event
        ###file name format: gw_latency_dict_S240501dr-2-Preliminary_test.pickle
        first_split = i.split('_')
        second_split = first_split[3].split('-')
        loading_dict["Name"] = second_split[0]

        ###collect timestamps of creation for mock files
        timestamp = os.path.getmtime(i)
        datestamp = datetime.fromtimestamp(timestamp, timezone.utc)
        loading_dict["time_stamp"] = datestamp

        ###is it the first or last file sent? not all events have 3 files,
        ###so determinging "last" will require more work in another function below
        if '-1-' in str(i):
            loading_dict["first"] = True
            loading_dict["last"] = False
        else:
            loading_dict["first"] = False
        if '-3-' in str(i):
            loading_dict["last"] = True
        if '-2-' in str(i):
            if len(glob.glob(path + f'PickledMocks/*{second_split[0]}-3-*'.format())) > 0:
                loading_dict["last"] = False
            else:
                loading_dict["last"] = True
        ###500 seconds of data must be collected before i3 can run, anything more is delay
        delay = 500 - mocks["Ligo_Latency"]
        loading_dict["we_had_to_wait"] = delay
        event_dict = event_dict.append(loading_dict, ignore_index=True)
    return event_dict
event_dict = sort_mocks(mock_files=mock_files)

###now to make the plots###
ed = event_dict

###call function incase mocks get interrupted
if (Time(datetime.now(timezone.utc)).mjd - max(Time(ed["time_stamp"]).mjd)) > 3600.*3/86400:
        dial_up()
        x = (Time(datetime.now(timezone.utc)).mjd - max(Time(ed["time_stamp"]).mjd))*24
        print("It has been " +str(x) +" hours since last update to gw mocks.")

tl = {"latency": [], "time_stamp": []} #total latency
il = {"latency": [], "time_stamp": []} #icecube latency
ll = {"latency": [], "time_stamp": []} #lvk latency
ts_list = []
###find the latencies associated with each institute, based only on the first map sent
for i in range(len(ed)):
    if ed["first"][i] == True:
        tl["latency"].append(ed["Total_Latency"][i])
        tl["time_stamp"].append(ed["time_stamp"][i])
        il["latency"].append(ed["IceCube_Latency"][i])
        il["time_stamp"].append(ed["time_stamp"][i])
        ll["latency"].append(ed["LVK_Latency"][i])
        ll["time_stamp"].append(ed["time_stamp"][i])
        ###find the events that have background trials run on them (increases latency)
        def find_ts(i, ts_list):
            ev_name = ed["Name"][i]
            find = glob.glob(out_put+f'/202*{ev_name}*/TS_distribution.png')
            if len(find)>0:
                ts_list.append(True)
            else:
                ts_list.append(False)
            return ts_list
        ts_list = find_ts(i=i, ts_list=ts_list)
    else:
        continue
ts_list = np.array(ts_list)

###tally amount of reports per day
unique_days = []
total = []
t = 0
preliminary = []
p = 0
initial = []
init = 0

for k in range(1, len(ed["time_stamp"])):
    if ed["time_stamp"][k-1].day == ed["time_stamp"][k].day:
        id = ed["Name"][k]
        if ed["Trigger_Time"][k-1] == ed["Trigger_Time"][k]:
            continue
        else:
            if len(glob.glob(out_put + f'/202*{id}-1-Preliminary_test/GFU_rate_plot.png')) > 0:
                p += 1
                t += 1
            if len(glob.glob(out_put + f'/202*{id}-2-Preliminary_test/GFU_rate_plot.png')) > 0:
                p += 1
                t += 1
            if len(glob.glob(out_put + f'/202*{id}-3-Initial_test/GFU_rate_plot.png')) > 0:
                init += 1
                t += 1
    else:
        unique_days.append(ed["time_stamp"][k-1].date())
        total.append(t)
        t = 0
        preliminary.append(p)
        p = 0
        initial.append(init)
        init = 0
print("Creating plots now")
fig, ax= plt.subplots(figsize=(12,6))

plt.title("Reports Per Day")
plt.xlabel("Date")
plt.ylabel("Number of Reports")

ax.plot_date(unique_days, total, color = 'black')
ax.plot_date(unique_days, preliminary, color = 'green')
ax.plot_date(unique_days, initial, color = 'orange')

now = datetime.now(timezone.utc)
past = now + relativedelta(months=-2)

ax.set_xlim(past, now)
ax.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
plt.xticks(rotation=18, ha="right")
plt.legend(["Total", "Preliminary", "Initial"], loc = "upper left")

save_path='/home/mromfoe/public_html/O4_followup_monitoring/ReportsPerDay_liveupdate.png'
plt.tight_layout()
fig.savefig(save_path)

###latency for the first report###
fig, ax= plt.subplots(figsize=(12,6))

plt.title("Latency per Report by Date")
plt.xlabel("Date")
plt.ylabel("Latency in Seconds")

ax.plot_date(tl["time_stamp"], tl["latency"], color = 'black')
ax.plot_date(il["time_stamp"], il["latency"], color = 'cyan')
ax.plot_date(ll["time_stamp"], ll["latency"], color = 'red')

now = datetime.now(timezone.utc)
past = now + relativedelta(weeks=-1)

ax.set_xlim(past, now)
ax.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
plt.xticks(rotation=18, ha="right")
plt.legend(["Total Latency", "IceCube Latency", "LVK Latency"], loc = "upper left")

save_path='/home/mromfoe/public_html/O4_followup_monitoring/Latencies_by_Date_liveupdate.png'
plt.tight_layout()
fig.savefig(save_path)

###Latencies by institute
n_bins = 25

fig, axs = plt.subplots(figsize = (10, 7))
plt.xlabel("Latency in Seconds: Last Bin is Overflow")
plt.ylabel("# of GCN Alerts")
plt.title("Total Latency: First Map")

total_med = round(np.median(np.array(tl["latency"])), 2)
bg_med = round(np.median(np.array(tl["latency"])[ts_list]), 2)
nbg_med = round(np.median(np.array(tl["latency"])[~ts_list]), 2)
tl["latency"] = np.array(tl["latency"])
tl["latency"][tl["latency"] > 3000] =  3000 ###creating overflow bin
h,b,_= axs.hist(tl["latency"], color = "c", bins = n_bins, label = "No BG Trial")
axs.hist((np.array(tl["latency"])[ts_list]), color = "orange", bins=b, label = "BG Trial")

percent_ts = []
pt = 0

for value in range(len(ts_list)):
        if ts_list[value] > 0:
                pt+=1
                percent_ts.append(pt)
percentage = round(len(percent_ts)/len(ts_list)*100, 2)

plt.axvline(total_med, color = 'r', linestyle = "dashdot", label = "Median: {:.2f}".format(total_med))
plt.axvline(bg_med, color = 'k', linestyle = 'dotted', label = "BG Trial Median: {:.2f}".format(bg_med))
plt.axvline(nbg_med, color = 'k', linestyle = 'dashed', label = "No BG Trial Median: {:.2f}".format(nbg_med))
plt.axvline(1, 0, linestyle = "none", label = "Percent w/ BG Trials: {:.2f} %".format(percentage))
plt.legend(fontsize = '15')
axs.grid(True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

save_path='/home/mromfoe/public_html/O4_followup_monitoring/Total_Latency_liveupdate.png'
plt.tight_layout()
fig.savefig(save_path)

###i3
fig, axs = plt.subplots(figsize = (10, 7))
plt.xlabel("Latency in Seconds: Last Bin is Overflow")
plt.ylabel("# of GCN Alerts")
plt.title("Latency on IceCube's Side: First Map")

total_med = round(np.median(np.array(il["latency"])), 2)
bg_med = round(np.median(np.array(il["latency"])[ts_list]), 2)
nbg_med = round(np.median(np.array(il["latency"])[~ts_list]), 2)
il["latency"] = np.array(il["latency"])
il["latency"][il["latency"] > 1200] =  1200 ###creating overflow bin
h,b,_= axs.hist(il["latency"], color = "c", bins = n_bins, label = "No BG Trial")
axs.hist((il["latency"])[ts_list], color = "orange", bins=b, label = "BG Trial")

plt.axvline(total_med, color = 'r', linestyle = "dashdot", label = "Median: {:.2f}".format(total_med))
plt.axvline(bg_med, color = 'k', linestyle = 'dotted', label = "BG Trial Median: {:.2f}".format(bg_med))
plt.axvline(nbg_med, color = 'k', linestyle = 'dashed', label = "No BG Trial Median: {:.2f}".format(nbg_med))
plt.axvline(1, 0, linestyle = "none", label = "Percent w/ BG Trials: {:.2f} %".format(percentage))
plt.legend(fontsize = '15')
axs.grid(True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

save_path='/home/mromfoe/public_html/O4_followup_monitoring/First_IceCube_Latency_liveupdate.png'
plt.tight_layout()
fig.savefig(save_path)

###lvk
fig, axs = plt.subplots(figsize = (10, 7))
plt.xlabel("Latency in Seconds: Last Bin is Overflow")
plt.ylabel("# of GCN Alerts")
plt.title("Latency on LVK's Side: First Map")

total_med = round(np.median(np.array(ll["latency"])), 2)
ll["latency"] = np.array(ll["latency"])
ll["latency"][ll["latency"] > 1200] =  1200 ###creating overflow bin
h,b,_= axs.hist(ll["latency"], color = "r", bins = n_bins)

plt.axvline(total_med, color = 'k', linestyle = "dashdot", label = "Median: {:.2f}".format(total_med))
plt.legend(fontsize = '15')
axs.grid(True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

save_path='/home/mromfoe/public_html/O4_followup_monitoring/First_LVK_Latency_liveupdate.png'
plt.tight_layout()
fig.savefig(save_path)

###find percentage of events we had to wait for
wait=[]
w=0
first_delay = []
for f in range(len(ed["first"])):
    if ed["first"][f]==True:
        first_delay.append(ed["we_had_to_wait"][f])
        if ed["we_had_to_wait"][f] > 0:
                w+=1
                wait.append(w)
percentage = round(len(wait)/len(first_delay)*100, 2)
fig, axs = plt.subplots(figsize = (10, 7))
plt.xlabel("Seconds of Delay")
plt.ylabel("# of GCN Alerts")
plt.title("500 Second Delay - First Map")
h,b,_= axs.hist(first_delay, bins = np.linspace(-2000, 500, num=25))

plt.axvline(0, color = 'k', linestyle = "solid", label = "We waited for items to the righ")

plt.legend(['We waited for: {:.2f}'.format(percentage)+'%'], fontsize = '18')

axs.grid(True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

save_path='/home/mromfoe/public_html/O4_followup_monitoring/500_sec_delay_liveupdate.png'

plt.tight_layout()
fig.savefig(save_path)

###text files for webpage
fig = plt.figure(figsize=(3.70, .20))
ax = fig.add_axes([0, 0, 1, 1])
plt.plot([0,5],[0,100],color='white')
ax.text(-0.1, 15, "Page Last Updated: {} UTC".format(now))
plt.savefig(save_path)

print("Creating tables")
df = pd.DataFrame({"Name": ed["Name"][-15::-1],
                    "Merger Time": ed["Trigger_Time"][-15::-1],
                    "GCN Alert": ed["GCN_Alert"][-15::-1],
                    "Script Finishes": ed["End_Time"][-15::-1],
                    "Total Latency in Seconds": ed["Total_Latency"][-15::-1]})
html1 = df.to_html()

text_file1 = open("/home/mromfoe/public_html/O4_followup_monitoring/Recent_Runs.html", "w")
text_file1.write(html1)
text_file1.close()
df = pd.DataFrame({"Name": [],
                "Merger Time": [],
                "GCN Alert": [],
                "Script Finishes": [],
                "Total Latency in Seconds": []})
for event in range(len(ed)):
    if ed["first"][event] == True:
        entry = {"Name": ed["Name"][event],
                "Merger Time": ed["Trigger_Time"][event],
                "GCN Alert": ed["GCN_Alert"][event],
                "Script Finishes": ed["End_Time"][event],
                "Total Latency in Seconds": ed["Total_Latency"][event]}
        df.loc[len(df)] = entry
    else:
        continue

html2 = df.to_html()

text_file2 = open("/home/mromfoe/public_html/O4_followup_monitoring/archive.html", "w")
text_file2.write(html2)
text_file2.close()

###make p-value distribution
print("making p-value plot")
def make_bg_pval_dist(fontsize=15, lower_y_bound=-3.5, load_all=False):
    # function to make pval dist. lower_y_bound arg gives the exponent to set the lower y-axis 
    # limit, e.g. 10^-3
    all_maps_saved_pkl=sorted(glob.glob('/data/user/jthwaites/o4-mocks/*/*.pickle'), reverse=True)
    saved_mock_pkl=[all_maps_saved_pkl[0]]

    for mock in all_maps_saved_pkl:
        event_name=mock.split('-')[-3]
        event_name=event_name.split('/')[-1]
        if event_name not in saved_mock_pkl[-1]:
            saved_mock_pkl.append(mock)

    all_mocks={}
    #if more than 1000 found - load most recent only (otherwise takes too long)  
    if len(saved_mock_pkl)>1000 and not load_all:
        print('Loading most recent 1000 mocks (may take a while)')
        saved_mock_pkl = saved_mock_pkl[0:1000]
    else:
        print('Loading %i mocks (may take a while)'%(len(saved_mock_pkl)))

    for mock in saved_mock_pkl:
        with open(mock,'rb') as f:
            result=pickle.load(f)
            all_mocks[result['name']]=result['p']
    print('Done loading mocks.')

    matplotlib.rcParams.update({'font.size':fontsize})
    plt.figure(figsize = (10,6), dpi=300)
    #ax.tick_params(labelsize=fontsize)

    p_x_vals = np.logspace(-4,0.,21)
    n, bins, patches = plt.hist([all_mocks[name] for name in all_mocks.keys()], 
                                weights = np.ones(len(all_mocks)) / len(all_mocks), bins = p_x_vals)
    
    lt_10per = sum(n[bins[:-1]<=0.1])
    lt_1per=sum(n[bins[:-1]<=0.01])
    
    #uniform_bins=np.logspace(lower_y_bound,0.,int(abs(lower_y_bound*7))+1) #evenly spaced bins in logspace
    #plt.step(uniform_bins[1:], np.diff(uniform_bins), label = 'Uniform p-value expectation', lw = 3.)
    plt.step(p_x_vals[1:], np.diff(p_x_vals), label = 'Uniform p-value distribution', lw = 3.)
    plt.plot([0.1,0.1], [10**lower_y_bound, 1e0],linestyle='dotted', label=f'{lt_10per*100.:.2f} \% of p-values $<$ 0.1')
    plt.plot([0.01, 0.01], [10**lower_y_bound, 1e0], linestyle='dashed',label=f'{lt_1per*100.:.2f} \% of p-values $<$ 0.01')

    plt.xscale('log')
    plt.yscale('log')
    plt.gca().invert_xaxis()
    plt.grid(which = 'both', alpha = 0.2)
    plt.xlim(1.1e0,1e-4)
    plt.ylim(10**lower_y_bound, 1e0)
    plt.xlabel('p-value', fontsize = fontsize)
    plt.ylabel('Fraction of Analyses', fontsize = fontsize)
    plt.legend(loc = 1, fontsize = fontsize)
    plt.title('{} Mock p-values as of {}'.format(len(saved_mock_pkl), str(date.today()),fontsize=20))

    save_path='/home/mromfoe/public_html/O4_followup_monitoring/mock_pvals_liveupdate.png'
    plt.savefig(save_path, dpi=200, bbox_inches='tight')
    print('Figure saved to file: ', save_path)
    
make_bg_pval_dist()