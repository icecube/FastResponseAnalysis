#!/usr/bin/env python
###For text files of timestamps & latencies from gw pipeline, go to end of script and call function
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams['axes.titlesize'] = 18
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16
import matplotlib.pyplot as plt
import numpy as np
import pickle, glob, os, datetime
from datetime import date
from statistics import median
#from statistics import mode
from matplotlib.dates import DateFormatter
from astropy.time import Time
import pandas as pd
#from IPython.display import HTML
#from IPython.display import display
import subprocess

def dial_up(who="jessie"):
    cell_tower = "/cvmfs/icecube.opensciencegrid.org/users/jthwaites/"
    #halp = "https://icecube.wisc.edu/~jthwaites/FastResponse/error_call.xml"
    subprocess.call([cell_tower+"make_call.py", f"--{who}=True", '--troubleshoot=True'])#, "--call_file", halp])
    #print([cell_tower+"make_call.py", f"--{who}=True", '--troubleshoot=True'])

#path = "/data/user/mromfoe/software/fastresponse/output/"
path = "/data/user/jthwaites/FastResponseAnalysis/output/"

out_put = '/data/user/jthwaites/o4-mocks/'

#Creating readable (and callable) files from ALL pickle files previously created in gw_gcn_listener
Pickle_to_text = sorted(glob.glob(path+'PickledMocks/*MS*.pickle'))

all_dictionary = {'Trigger_Time': [], 'GCN_Alert': [], 'End_Time': [],
                        'LVK_Latency': [], 'IceCube_Latency': [], 'Total_Latency': [],
                            'We_had_to_wait:': [], 'Name': [], 'Time_Stamp': []}

for file in Pickle_to_text:
        if os.path.exists(file):
                with open(file, 'rb') as pickle_now:
                        Entries = pickle.load(pickle_now)
                if int(Entries["Trigger_Time"]) == 59957:
                        continue #This is January 13th 2023, dates were messed up.
                for Key in Entries.keys():
                        if Key=="Ligo_Latency":
                                all_dictionary['LVK_Latency'].append(Entries[Key])
                        else:
                                all_dictionary[Key].append(Entries[Key])
                all_dictionary['Name'].append(file)
for file in all_dictionary["Name"]:
        date_format = "%Y-%m-%d"
        c_timestamp = os.path.getctime(file)
        c_datestamp = datetime.datetime.utcfromtimestamp(c_timestamp)
        all_dictionary['Time_Stamp'].append(c_datestamp)
print('Finished loading all latencies.')

#Setting conditions for call to be made if script fails. Condition: if more than 2 hours pass between alerts
if (Time(datetime.datetime.utcnow()).mjd - max(Time(all_dictionary["Time_Stamp"]).mjd)) > 3600.*2/86400:
        #print(datetime.datetime.utcnow())
        #print(max(all_dictionary["Time_Stamp"]))
        dial_up()
        x = (Time(datetime.datetime.utcnow()).mjd - max(Time(all_dictionary["Time_Stamp"]).mjd))*24
        print("It has been " +str(x) +" hours since last update to gw mocks.")
#        print("Uh oh... spaghetti-o's")

#Sorting pickle files by date created and by most correct version
Pickle_Files = sorted(glob.glob(path+'PickledMocks/*.pickle'))

date_format = "%Y-%m-%d"

Quality_Pickles = []

for file in Pickle_Files:
        c_timestamp = os.path.getctime(file)
        c_datestamp = datetime.datetime.fromtimestamp(c_timestamp)
        if (c_datestamp >= datetime.datetime.strptime('2022-12-05', date_format)):
                Quality_Pickles.append(file)

#Collecting the first maps created for each event for analysis
First_Batch={'Trigger_Time': [], 'GCN_Alert': [], 'End_Time': [],
                'LVK_Latency': [], 'IceCube_Latency': [], 'Total_Latency': [],
                'We_had_to_wait:': [], 'Name': [], 'Time_Stamp': []}        
First_Runs= sorted(glob.glob((path+'PickledMocks/*MS*-1-*.pickle')))

for file in First_Runs:
        if os.path.exists(file):
                with open(file, 'rb') as pickle_now:
                        Entries = pickle.load(pickle_now)
                        if int(Entries["Trigger_Time"]) == 59957:
                                continue
                        for Key in Entries.keys():
                                if Key=="Ligo_Latency":
                                        First_Batch['LVK_Latency'].append(Entries[Key])
                                else:
                                        First_Batch[Key].append(Entries[Key])        
                First_Batch['Name'].append(file)
for file in First_Batch["Name"]:
        date_format = "%Y-%m-%d"
        c_timestamp = os.path.getctime(file)
        c_datestamp = datetime.datetime.fromtimestamp(c_timestamp)
        First_Batch['Time_Stamp'].append(c_datestamp)
print('Finished loading 1st map latencies.')

#Breakdown full name of files created so only necessary parts are saved (i.e., the code for the event and number of the map)
pieces = [string.split('-') for string in Quality_Pickles]

Puzzle_Box = {}

for prefix in np.unique([pieces[i][0] for i in range(len(pieces))]):
        Puzzle_Box[prefix]=[[],[]]
for bit in pieces:
        Puzzle_Box[bit[0]][0].append(int(bit[1]))
        Puzzle_Box[bit[0]][1].append(bit[2])

answer = []
names = []
#Finds the most recent map created for each event. This is usually '2' or '3'. 
for key in Puzzle_Box.keys():
        high_file=max(Puzzle_Box[key][0])
        for i in range(len(Puzzle_Box[key][0])):
                if Puzzle_Box[key][0][i]==high_file:
                        answer.append(key+'-'+str(high_file)+'-'+Puzzle_Box[key][1][i])
                        names.append(key[-9:]+'-'+str(high_file)+'-'+Puzzle_Box[key][1][i].split('.')[0])

Last_Batch = {'Trigger_Time': [], 'GCN_Alert': [], 'End_Time': [],
                        'LVK_Latency': [], 'IceCube_Latency': [], 'Total_Latency': [],
                        'We_had_to_wait:': [], 'Name': [], 'Time_Stamp': []}

#Make keys in dictionary callable
for file in answer:
    if os.path.exists(file):
        with open(file, 'rb') as pickle_in:
            Brine = pickle.load(pickle_in)
            if int(Brine["Trigger_Time"]) == 59957: #This is January 13, 2023.
                    continue                        #There was an error causing multiple files to be falsely saved to this date.
            for Key in Brine.keys():
                if Key=="Ligo_Latency":
                        Last_Batch['LVK_Latency'].append(Brine[Key])
                else:
                        Last_Batch[Key].append(Brine[Key])
        Last_Batch['Name'].append(file)
for file in Last_Batch["Name"]:
        date_format = "%Y-%m-%d"
        c_timestamp = os.path.getctime(file)
        c_datestamp = datetime.datetime.fromtimestamp(c_timestamp)
        Last_Batch['Time_Stamp'].append(c_datestamp)
print('Finished loading last map latencies')

First_TS = []
Last_TS = []

def find_ts(Batch, TS_List):

        split_names = [string.split('/') for string in Batch]

        split_bits = []
        prime_split = []
        final_bit = []

        for i in range(len(split_names)):
                split_bits.append(split_names[i][-1])

        middle_split = [string.split('.') for string in split_bits]

        for i in range(len(middle_split)):
                prime_split.append(middle_split[i][-2])

        chop = [string.split('_') for string in prime_split]

        for i in range(len(chop)):
                final_bit.append(chop[i][-2]+'_'+chop[i][-1])

        for i in final_bit:
                FTS=glob.glob(out_put+f'/202*{i}/TS_distribution.png')
                if len(FTS)>0:
                        TS_List.append(True)
                else:
                        TS_List.append(False)
        TS_List = np.array(TS_List)
        return TS_List
First_TS = find_ts(Batch = First_Batch["Name"], TS_List = First_TS)
Last_TS = find_ts(Batch = Last_Batch["Name"], TS_List = Last_TS)

#Find the outliers (those that take a long time to run/process)
outliers_ice = []
outliers_tot = []

for item in range(len(First_Batch["IceCube_Latency"])):
        if First_Batch["IceCube_Latency"][item]>=1200:
                outliers_ice.append(First_Batch["Name"][item])

for item in range(len(First_Batch["Total_Latency"])):
        if First_Batch["Total_Latency"][item]>=3000:
                outliers_tot.append(First_Batch["Name"][item])

icecube_outliers = []
total_outliers = []

def outlier_name(outlier_list, outlier_names):
        nb = [string.split('/') for string in outlier_list]

        out_bit=[]

        for i in range(len(nb)):
                out_bit.append(nb[i][-1])

        last_bit = [string.split('_') for string in out_bit]

        for i in range(len(last_bit)):
                outlier_names.append(last_bit[i][-2])
outlier_name(outlier_list = outliers_ice, outlier_names = icecube_outliers)
outlier_name(outlier_list = outliers_tot, outlier_names = total_outliers)

#creating dictionaries of outliers so they can be quickly and easily identified
d_outlier_ice = {"Name": [], "Trigger Time": [], "End Time": []}
d_outlier_tot = {"Name": [], "Trigger Time": [], "End Time": []}

for n in icecube_outliers:
        for i in range(len(First_Batch["Name"])):
                if n in First_Batch["Name"][i]:
                        d_outlier_ice["Name"].append(n)
                        d_outlier_ice["Trigger Time"].append(First_Batch["Trigger_Time"][i])
                        d_outlier_ice["End Time"].append(First_Batch["End_Time"][i])
for n in total_outliers:
        for i in range(len(First_Batch["Name"])):
                if n in First_Batch["Name"][i]:
                        d_outlier_tot["Name"].append(n)
                        d_outlier_tot["Trigger Time"].append(First_Batch["Trigger_Time"][i])
                        d_outlier_tot["End Time"].append(First_Batch["End_Time"][i])

a = all_dictionary
###Create pandas table to display on website
#This table shows only the most recently received map data
name_bits = [string.split('/') for string in all_dictionary['Name']]

new_names=[]
prime_bits=[]

for i in range(len(name_bits)):
        prime_bits.append(name_bits[i][-1])

last_bit = [string.split('_') for string in prime_bits]

for i in range(len(prime_bits)):
        new_names.append(last_bit[i][-2])

Merger=[]
for i in range(len(a["Trigger_Time"])):
        merg= Time(a["Trigger_Time"][i], format="mjd").iso
        Merger.append(merg)
Alert=[]
for i in range(len(a["GCN_Alert"])):
        aler= Time(a["GCN_Alert"][i], format="mjd").iso
        Alert.append(aler)
End=[]
for i in range(len(a["End_Time"])):
        en= Time(a["End_Time"][i], format="mjd").iso
        End.append(en)

total_late=[]
for i in range(len(a["Total_Latency"])):
        tot= a["Total_Latency"][i]%3600
        total_late.append(tot)

n_names= new_names[:-16:-1]
m_merger= Merger[:-16:-1]
a_alert= Alert[:-16:-1]
e_end= End[:-16:-1]
tot_late= total_late[:-16:-1]

df = pd.DataFrame({"Name": n_names,
                        "Merger Time": m_merger,
                        "GCN Alert": a_alert,
                        "Script Finishes": e_end,
                        "Total Latency in Seconds": tot_late})
html1 = df.to_html()

text_file1 = open("/home/mromfoe/public_html/O4_followup_monitoring/Recent_Runs.html", "w")
text_file1.write(html1)
text_file1.close()

###This table displays the LAST map from each event
name_bits = [string.split('/') for string in First_Batch['Name']]

FB=First_Batch

new_names=[]
prime_bits=[]

for i in range(len(name_bits)):
        prime_bits.append(name_bits[i][-1])

nlast_bit = [string.split('_') for string in prime_bits]

for i in range(len(prime_bits)):
        new_names.append(nlast_bit[i][-2])

Merger=[]
for i in range(len(FB["Trigger_Time"])):
        merg= Time(FB["Trigger_Time"][i], format="mjd").iso
        Merger.append(merg)
Alert=[]
for i in range(len(FB["GCN_Alert"])):
        aler= Time(FB["GCN_Alert"][i], format="mjd").iso
        Alert.append(aler)
End=[]
for i in range(len(FB["End_Time"])):
        en= Time(FB["End_Time"][i], format="mjd").iso
        End.append(en)

total_late=[]
for i in range(len(FB["Total_Latency"])):
        tot= FB["Total_Latency"][i]
        total_late.append(tot)

df = pd.DataFrame({"Name": new_names,
                        "Merger Time": Merger,
                        "GCN Alert": Alert,
                        "Script Finishes": End,
                        "Total Latency in Seconds": total_late})

html2 = df.to_html()

text_file2 = open("/home/mromfoe/public_html/O4_followup_monitoring/archive.html", "w")
text_file2.write(html2)
text_file2.close()

#Create histograms of time data collected
n_bins = 25
###
fig, axs = plt.subplots(figsize = (10, 7))
plt.xlabel("Latency in Seconds: Last Bin is Overflow")
plt.ylabel("# of GCN Alerts")
plt.title("Latency on LVK's Side: First Map")
med_LVK = round(median(np.array(First_Batch['LVK_Latency'])), 2)
First_Batch["LVK_Latency"] = np.array(First_Batch["LVK_Latency"])
First_Batch["LVK_Latency"][First_Batch["LVK_Latency"] > 1200] =  1200
h,b,_= axs.hist(First_Batch['LVK_Latency'], color = "r", bins = n_bins)

tallest = max(h)

plt.axvline(med_LVK, color = 'k', linestyle = 'dashed', label = "LVK Median")
plt.text(med_LVK*1.3, tallest*0.8, 'Median: {:.2f}'.format(med_LVK), fontsize = "15")

plt.legend(fontsize = "18")

axs.grid(b = True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

save_path='/home/mromfoe/public_html/O4_followup_monitoring/First_LVK_Latency_liveupdate.png'
plt.tight_layout()
fig.savefig(save_path)
fig.savefig("LVK_Latency.png")
plt.close()

from PIL import Image, ImageDraw, ImageFont

img = Image.new('RGB', (370, 20), "white")

d1 = ImageDraw.Draw(img)

#fontname = "/data/user/mromfoe/software/fastresponse/fast_response/MonitoringAndMocks/A101HLVN.ttf"
fontname = "/home/mromfoe/public_html/O4_followup_monitoring/A101HLVN.ttf"
fontsize = 16

MyFont = ImageFont.truetype(fontname, fontsize)

now = Time(datetime.datetime.utcnow()).iso

d1.text((2, 0), "Page Last Updated: {} UTC".format(now), 
        fill = (0, 0, 0), font = MyFont)

save_path='/home/mromfoe/public_html/O4_followup_monitoring/Update_Time.png'

img.save(save_path)
img.save("Update_Time.png")
###
fig, axs = plt.subplots(figsize = (10, 7))
plt.xlabel("Latency in Seconds")
plt.ylabel("# of GCN Alerts")
plt.title("Latency on IceCube's Side: First Map")
med_red = round(median(np.array(First_Batch['IceCube_Latency'])), 2)
med_blue = round(median(np.array(First_Batch['IceCube_Latency'])[~First_TS]), 2)
med_orange = round(median(np.array(First_Batch['IceCube_Latency'])[First_TS]), 2)
First_Batch["IceCube_Latency"] = np.array(First_Batch["IceCube_Latency"])
First_Batch["IceCube_Latency"][First_Batch["IceCube_Latency"] > 1200] =  1200
h,b,_= axs.hist(First_Batch['IceCube_Latency'], color = "c", bins = n_bins, label = "No BG Trial")
axs.hist(np.array(First_Batch['IceCube_Latency'])[First_TS], color = "orange", bins = b, label = "BG Trial")

percent_ts = []
pt = 0

for value in range(len(First_TS)):
        if First_TS[value] > 0:
                pt+=1
                percent_ts.append(pt)
percentage = round(len(percent_ts)/len(First_TS)*100, 2)

plt.axvline(med_red, color = 'r', linestyle = "dashdot", label = "Median")
plt.axvline(med_blue, color = 'k', linestyle = 'dashed', label = "No BG Trial Median")
plt.axvline(med_orange, color = 'k', linestyle = 'dotted', label = "BG Trial Median")
plt.axvline(1, 0, linestyle = "none", label = "Percent w/ BG Trials: {:.2f} %".format(percentage))

tallest = max(h)

plt.text(med_blue*0.15, tallest*0.8, 'No BG Median: {:.2f}'.format(med_blue), fontsize = "13")
plt.text(med_orange*1.05, tallest*0.5, 'BG Median: {:.2f}'.format(med_orange), fontsize = "13")
plt.text(med_red*1.1, tallest*0.55, 'Median: {:.2f}'.format(med_red), fontsize = "13")


plt.legend(fontsize = '15')

axs.grid(b = True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

save_path='/home/mromfoe/public_html/O4_followup_monitoring/First_IceCube_Latency_liveupdate.png'

plt.tight_layout()
fig.savefig(save_path)
fig.savefig("First_IceCube_Latency.png")
fig, axs = plt.subplots(figsize = (10, 7))
plt.close()
###
fig, axs = plt.subplots(figsize = (10, 7))
plt.xlabel("Latency in Seconds")
plt.ylabel("# of GCN Alerts")
plt.title("Latency on IceCube's Side: Last Map")
med_red = round(median(np.array(Last_Batch['IceCube_Latency'])), 2)
med_blue = round(median(np.array(Last_Batch['IceCube_Latency'])[~Last_TS]), 2)
med_orange = round(median(np.array(Last_Batch['IceCube_Latency'])[Last_TS]), 2)
Last_Batch["IceCube_Latency"] = np.array(Last_Batch["IceCube_Latency"])
Last_Batch["IceCube_Latency"][Last_Batch["IceCube_Latency"] > 600] =  600
h,b,_= axs.hist(Last_Batch['IceCube_Latency'], color = "c", bins= n_bins, label = "No BG Trial")
axs.hist(np.array(Last_Batch['IceCube_Latency'])[Last_TS], color = "orange", bins= b, label = "BG Trial")

percent_ts = []
pt = 0

for value in range(len(First_TS)):
        if First_TS[value] > 0:
                pt+=1
                percent_ts.append(pt)
percentage = round(len(percent_ts)/len(First_TS)*100, 2)

plt.axvline(med_red, color = 'r', linestyle = "dashdot", label = "Median")
plt.axvline(med_blue, color = 'k', linestyle = 'dashed', label = "No BG Trial Median")
plt.axvline(med_orange, color = 'k', linestyle = 'dotted', label = "BG Trial Median")
plt.axvline(1, 0, linestyle = "none", label = "Percent w/ BG Trials: {:.2f} %".format(percentage))

tallest = max(h)

plt.text(-20, tallest*0.8, 'No BG Median: {:.2f}'.format(med_blue), fontsize = "13")
plt.text(med_orange*1.05, tallest*0.5, 'BG Median: {:.2f}'.format(med_orange), fontsize = "13")
#plt.text(400, 150, "# of Background Trials: {:.2f}".format(sum(Last_TS)), fontsize = "15")
plt.text(med_red*1.1, tallest*0.7, 'Median: {:.2f}'.format(med_red), fontsize = "13")

plt.legend(fontsize = '18')

axs.grid(b = True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

save_path='/home/mromfoe/public_html/O4_followup_monitoring/Last_IceCube_Latency_liveupdate.png'

plt.tight_layout()
fig.savefig(save_path)
fig.savefig("Last_IceCube_Latency.png")
plt.close()
###

fig, axs = plt.subplots(figsize = (10, 7),
                        tight_layout = True)

plt.xlabel("Latency in Seconds")
plt.ylabel("# of GCN Alerts")
plt.title("Total Latency: First Map")
med_red = round(median(np.array(First_Batch['Total_Latency'])), 2)
med_blue = round(median(np.array(First_Batch['Total_Latency'])[~First_TS]), 2)
med_orange = round(median(np.array(First_Batch['Total_Latency'])[First_TS]), 2)
First_Batch["Total_Latency"] = np.array(First_Batch["Total_Latency"])
First_Batch["Total_Latency"][First_Batch["Total_Latency"] > 3000] =  3000
h,b,_= axs.hist(First_Batch['Total_Latency'], color = "c", bins = n_bins, label = "No BG Trial")
axs.hist(np.array(First_Batch['Total_Latency'])[First_TS], color = "orange", bins=b, label = "BG Trial")

percent_ts = []
pt = 0

for value in range(len(First_TS)):
        if First_TS[value] > 0:
                pt+=1
                percent_ts.append(pt)
percentage = round(len(percent_ts)/len(First_TS)*100, 2)

plt.axvline(med_red, color = 'r', linestyle = "dashdot", label = "Median")
plt.axvline(med_blue, color = 'k', linestyle = 'dashed', label = "No BG Trial Median")
plt.axvline(med_orange, color = 'k', linestyle = 'dotted', label = "BG Trial Median")
plt.axvline(1, 0, linestyle = "none", label = "Percent w/ BG Trials: {:.2f} %".format(percentage))

tallest = max(h)

plt.text(-100, tallest*0.8, 'No BG Med.: {:.2f}'.format(med_blue), fontsize = "13")
plt.text(med_orange*1.1, tallest*0.6, 'BG Med.: {:.2f}'.format(med_orange), fontsize = "13")
plt.text(med_red*1.1, tallest*0.7, 'Median: {:.2f}'.format(med_red), fontsize = "13")
#plt.text(4000, 175, "# of BG Trials: {:.2f}".format(sum(First_TS)), fontsize = "15")

plt.legend(fontsize = '15')

save_path='/home/mromfoe/public_html/O4_followup_monitoring/Total_Latency_liveupdate.png'

axs.grid(b = True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

plt.tight_layout()
fig.savefig(save_path)
fig.savefig("First_Total_Latency.png")

###
wait=[]
w=0
for value in range(len(First_Batch["We_had_to_wait:"])):
        if First_Batch["We_had_to_wait:"][value] > 0:
                w+=1
                wait.append(w)
percentage = round(len(wait)/len(First_Batch["We_had_to_wait:"])*100, 2)

fig, axs = plt.subplots(figsize = (10, 7))
plt.xlabel("Seconds of Delay")
plt.ylabel("# of GCN Alerts")
plt.title("500 Second Delay")
h,b,_= axs.hist(First_Batch['We_had_to_wait:'], bins = np.linspace(-2000, 500, num=25))

plt.axvline(0, color = 'k', linestyle = "solid", label = "We waited for items to the righ")

plt.legend(['We waited for: {:.2f}'.format(percentage)+'%'], fontsize = '18')

axs.grid(b = True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

save_path='/home/mromfoe/public_html/O4_followup_monitoring/500_sec_delay_liveupdate.png'

plt.tight_layout()
fig.savefig(save_path)
fig.savefig("500_sec_delay.png")

###Comparison of latencies for each day.
fig, ax= plt.subplots(figsize=(12,6))

ax.plot_date(First_Batch["Time_Stamp"], First_Batch["LVK_Latency"], color = 'red')
ax.plot_date(First_Batch["Time_Stamp"], First_Batch["IceCube_Latency"], color = 'cyan')
ax.plot_date(First_Batch["Time_Stamp"], First_Batch["Total_Latency"], color = 'black')

ax.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
plt.xticks(rotation=18, ha="right")

plt.ylim(0, 6000)

plt.title("Latencies by Date")
plt.legend(["LVK Latency", "IceCube Latency", "Total Latency"])
plt.xlabel("Date")
plt.ylabel("Latency in Seconds")

save_path='/home/mromfoe/public_html/O4_followup_monitoring/Latencies_by_Date_liveupdate.png'

plt.tight_layout()
fig.savefig(save_path)
fig.savefig("Latencies_by_Date.png")

###Zoomed in view

fig, ax= plt.subplots(figsize=(12,6))

ax.plot_date(First_Batch["Time_Stamp"], First_Batch["LVK_Latency"], color = 'red')
ax.plot_date(First_Batch["Time_Stamp"], First_Batch["IceCube_Latency"], color = 'cyan')
ax.plot_date(First_Batch["Time_Stamp"], First_Batch["Total_Latency"], color = 'black')

ax.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')
plt.xticks(rotation=18, ha="right")

plt.ylim(0, 2000)

plt.title("Latencies by Date (Zoomed in)")
plt.legend(["LVK Latency", "IceCube Latency", "Total Latency"])
plt.xlabel("Date")
plt.ylabel("Latency in Seconds")

save_path='/home/mromfoe/public_html/O4_followup_monitoring/Zoomed_Latencies_by_Date_liveupdate.png'

plt.tight_layout()
fig.savefig(save_path)
fig.savefig("Zoomed_Latencies_by_Date.png")

###Number of reports per day
a = all_dictionary

perday = []
unique_days = []
init = []
prelim = []
update = []
n = 0
t = 0
p = 0
u = 0
for i in range(1, len(all_dictionary["Time_Stamp"])):
        if a["Time_Stamp"][i-1].day == a["Time_Stamp"][i].day:
                n += 1
                if "Preliminary" in a["Name"][i]:
                        p += 1
                if "Initial" in a["Name"][i]:
                        t += 1
                if "Update" in a["Name"][i]:
                        u += 1
                
        else:
                unique_days.append(a["Time_Stamp"][i-1].date())
                perday.append(n)
                n = 0
                init.append(t)
                t = 0
                prelim.append(p)
                p = 0
                update.append(u)
                u = 0

fig, ax= plt.subplots(figsize=(12,6))

ax.plot_date(unique_days, perday, color = 'black')
ax.plot_date(unique_days, prelim, color = 'cyan')
ax.plot_date(unique_days, init, color = 'green')
ax.plot_date(unique_days, update, color = 'red')

ax.fmt_xdata = DateFormatter('%Y-%m-%d')
plt.xticks(rotation=18, ha="right")

plt.legend(["Total", "Preliminary", "Initial", "Update"], loc = "upper left")

plt.title("Reports per Day")
plt.xlabel("Date")
plt.ylabel("Number of Reports")

save_path='/home/mromfoe/public_html/O4_followup_monitoring/ReportsPerDay_liveupdate.png'

fig.savefig(save_path)
fig.savefig("ReportsPerDay.png")

###Find average time between mocks/runs

time_average = []
for i in range(len(First_Batch["Trigger_Time"]))[:-1]:
        delta_t = First_Batch["Trigger_Time"][i+1] - First_Batch["Trigger_Time"][i]
        time_average.append(delta_t*86400)

fig, axs = plt.subplots(figsize = (10, 7))

plt.xlabel("Time Between Runs in Seconds")
plt.ylabel("# of Runs")
plt.title("Average Delta T")

time_average = np.array(time_average)
time_average[time_average > 8000] =  8000
h,b,_= axs.hist(time_average, bins= np.linspace(0, 8000, 25))
med_delta = median(time_average)

tallest = max(h)

plt.axvline(med_delta, color = 'k', linestyle = 'dashed', label = "Median Delta T")
plt.text(med_delta*1.3, tallest*0.8, 'Median: {:.2f}'.format(med_delta), fontsize = "15")

plt.legend(fontsize = "18")

axs.grid(b = True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

plt.tight_layout()
fig.savefig("Linear_Delta_T.png")

time_average_too = []
for i in range(len(First_Batch["Time_Stamp"]))[:-1]:
        delta_t = Time(First_Batch["Time_Stamp"][i+1]).mjd - Time(First_Batch["Time_Stamp"][i]).mjd
        time_average_too.append(delta_t*86400)

fig, axs = plt.subplots(figsize = (10, 7))

plt.xlabel("Time Between Runs in Seconds")
plt.ylabel("# of Runs")
plt.title("Average Delta T")

time_average_too = np.array(time_average_too)
time_average_too[time_average_too > 8000] =  8000
h,b,_= axs.hist(time_average_too, bins= np.linspace(0, 8000, 25))
med_delta = median(time_average_too);

tallest = max(h)

plt.axvline(med_delta, color = 'k', linestyle = 'dashed', label = "Median Delta T")
plt.text(med_delta*1.3, tallest*0.8, 'Median: {:.2f}'.format(med_delta), fontsize = "15")

plt.legend(fontsize = "18")

axs.grid(b = True, color ='black',
        linestyle ='-.', linewidth = 0.5,
        alpha = 0.7)

plt.tight_layout()
fig.savefig("Linear_Delta_T_TIMESTAMP.png")
###p value plots
import matplotlib as mpl
mpl.use('agg')

def make_bg_pval_dist(fontsize=15, lower_y_bound=-3.5):
    # function to make pval dist. lower_y_bound arg gives the exponent to set the lower y-axis 
    # limit, e.g. 10^-3
    all_maps_saved_pkl=sorted(glob.glob('/data/user/jthwaites/o4-mocks/*/*.pickle'))[::-1]
    saved_mock_pkl=[all_maps_saved_pkl[0]]

    for mock in all_maps_saved_pkl:
        event_name=mock.split('-')[-3]
        event_name=event_name.split('/')[-1]
        if event_name not in saved_mock_pkl[-1]:
            saved_mock_pkl.append(mock)

    all_mocks={}
    print('Loading %i mocks (may take a while)'%(len(saved_mock_pkl)))
    for mock in saved_mock_pkl:
        with open(mock,'rb') as f:
            result=pickle.load(f)
            all_mocks[result['name']]=result['p']
    print('Done loading mocks.')

    mpl.rcParams.update({'font.size':fontsize})
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
    plt.savefig("mock_pvals.png")
    print('Figure saved to file: ', save_path)
    
make_bg_pval_dist()

def readtimes():
        a = all_dictionary

        file_object = open('MilestoneTimes.txt', "w+")
        file_object.write('\n' +"Trigger Time=" +repr(a["Trigger_Time"]) +'\n' +"GCN Alert=" +repr(a["GCN_Alert"]) 
                        +'\n'  +"End Time=" +repr(a["End_Time"]))
        file_object.close()

        file_object = open('GWLatency.txt', "w+")
        file_object.write('\n' +"LVK Latency=" +repr(a["LVK_Latency"]) +'\n' +"IceCube Latency=" 
                +repr(a["IceCube_Latency"]) +'\n'+ "Total Latency=" +repr(a["Total_Latency"]) +'\n' 
                +"We had to wait..." +repr(a["We_had_to_wait:"]) +"seconds." +'\n')
        file_object.close()
