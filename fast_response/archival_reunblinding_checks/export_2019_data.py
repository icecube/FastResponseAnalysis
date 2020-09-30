from skylab.datasets.GFUOnline_v001p00 import GFUOnline_v001p00
import numpy as np
import json
import icecube.realtime_gfu.eventcache
import icecube.realtime_tools.live
from icecube.phys_services import goodrunlist

def clean_json_events(evs):
    stream = icecube.realtime_gfu.eventcache.NeutrinoEventStream()
    events = stream.parse(evs)
    dtype = np.load('/data/ana/analyses/gfu_online/current/IC86_2018_data.npy').dtype
    eventlist = []
    for event in events:
        reordered = []
        for n in dtype.names:
            if not n in event.dtype.names: reordered.append(-1)
            else: reordered.append(event[n])
        eventlist.append(tuple(reordered))
    events = np.array(eventlist, dtype=dtype)
    return events

def format_db_runs(run_table, events, i3_grl):
    # loop through good runs to grab events & info
    i3_runs = i3_grl.keys()
    run_info = []
    exp = []
    for run in run_table:
        if run['stop'] is None:
            continue
        if run['run_number'] not in i3_runs:
            continue
        if i3_grl[run['run_number']]['active_strings'] < 80:
            continue
        if not run['latest_snapshot']['good_i3']:
            continue
        
        # start & stop time of run in mjd days
        run_start = run['start']
        run_stop = run['stop']

        mask = ((run_start <= events['time']) &
                (events['time'] <= run_stop) &
                (events['run'] == run['run_number']))

        exp.append(events[mask].copy())

        # append GRL info
        run_info.append([run['run_number'], run_start, run_stop,
                         exp[-1].size])

    # build the good run list
    dtype = [('run', int), ('start', float), ('stop', float),
             ('livetime', float), ('events', int), ('good_i3', bool)]
    grl = np.empty(len(run_info), dtype)
    run_info = np.transpose(run_info)
    grl['run'] = run_info[0]
    grl['start'] = run_info[1]
    grl['stop'] = run_info[2]
    grl['livetime'] = grl['stop'] - grl['start']
    grl['events'] = run_info[3]
    grl['good_i3'] = True

    # stick all good events in a single array
    exp = np.concatenate(exp)

    return exp, grl, grl['livetime'].sum()


if __name__=='__main__':
    official_grl = goodrunlist.GRL()

    with open('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/all_realtimeEvent.json', 'r') as f:
            docs = json.loads(f.read())
            
    exp_compare = np.load('/data/ana/analyses/gfu_online/current/IC86_2018_data.npy')
    grl_compare = np.load('/data/ana/analyses/gfu_online/current/GRL/IC86_2018_data.npy')

    min_compare_run = 131265 #not just np.min because I wanted to ignore test runs
    max_compare_run = np.max(grl_compare['run'])
    compare_docs = []
    new_docs = []
    for doc in docs:
        #Only keep runs from IC2018 to compare, throw out the early test runs
        if (doc['value']['data']['run_id'] >= min_compare_run) and (doc['value']['data']['run_id'] <= max_compare_run):
            compare_docs.append(doc)
        elif doc['value']['data']['run_id'] > max_compare_run:
            new_docs.append(doc)
        else:
            continue #There are a few events Michael gave from even before the IC2018 files      

    exp = clean_json_events(new_docs)
    start = np.min(exp['time']); stop = np.max(exp['time'])

    run_table = GFUOnline_v001p00.query_db_runs(start, stop)

    exp, grl, livetime = format_db_runs(run_table, exp, official_grl)

    exp.sort(order='time')
    grl.sort(order='run')
    
    trimmed_grl = grl[1:] #There will be one overlapping run at the beginning that's acutally in 2018

    np.save('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/2019_data/IC86_2019_data.npy', exp)
    np.save('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/2019_data/GRL/IC86_2019_data.npy', trimmed_grl)
