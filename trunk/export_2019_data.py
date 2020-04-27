from skylab.datasets.GFUOnline_v001p00 import GFUOnline_v001p00
import numpy as np
import json
import icecube.realtime_gfu.eventcache
import icecube.realtime_tools.live

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

def format_db_runs(run_table, events):
    #print(len(events))
    # loop through good runs to grab events & info
    run_info = []
    exp = []
    for run in run_table:

        # skip currently unfinished runs
        if run['stop'] is None:
            continue
        # skip bad runs
        if run['OK'] != 'OK':
            continue

        # start & stop time of run in mjd days
        run_start = run['start']
        run_stop = run['stop']

        mask = ((run_start <= events['time']) &
                (events['time'] <= run_stop) &
                (events['run'] == run['run_number']))
        print(len(events[mask]))
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

    # trim to ensure events and runs are between start/stop
    #exp = exp[(start <= exp['time']) & (exp['time'] <= stop)]
    #grl = grl[(start <= grl['stop']) & (grl['start'] <= stop)]
    #if grl['start'][0] < start:
    #    grl['start'][0] = start
    #    grl['livetime'][0] = grl['stop'][0] - grl['start'][0]
    #if grl['stop'][-1] > stop:
    #    grl['stop'][-1] = stop
    #    grl['livetime'][-1] = grl['stop'][-1] - grl['start'][-1]

    return exp, grl, grl['livetime'].sum()


if __name__=='__main__':
    with open('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/all_realtimeEvent.json', 'r') as f:
        docs = json.loads(f.read())
    exp = clean_json_events(docs[:1000])
    print(len(exp))
    print(np.unique(exp['run']))
    start = np.min(exp['time']); stop = np.max(exp['time'])

    run_table = GFUOnline_v001p00.query_db_runs(start, stop)    
    exp, grl, livetime = format_db_runs(run_table, exp)

    exp.sort(order='time')
    grl.sort(order='run')

    np.save('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/IC_2019_data.npy', exp)
    np.save('/data/user/apizzuto/fast_response_skylab/fast-response/trunk/IC_2019_grl.npy', grl)

    
