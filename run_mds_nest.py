import sys, pprint
import nest
import numpy as np
import matplotlib.pyplot as plt
import time
import ipdb
nest.set_verbosity('M_ALL')

def mds_nest_sim(stim, curr_arr, d_dt):
    time_stamps = np.arange(0,curr_arr.size*d_dt,d_dt)
    if stim == 'corrcost':
        sim_lenght = time_stamps[-1]
        n_time_steps = int(sim_lenght/d_dt)
        sim_type = 'corrcost'    

    elif stim == 'corrcostatratti':
        sim_lenght = time_stamps[-1]
        n_time_steps = int(sim_lenght/d_dt)
        sim_type = 'corrcostatratti'


    time_step = d_dt
    sim_time = time_step*n_time_steps
    t_arr = time_stamps 

    nest.ResetKernel()

    # Set simulation parameters
    nest.resolution = d_dt
    
    # create a step current generator
    cur_gen = nest.Create("step_current_generator")
    # set times at which current changes and amplitudes of step current
    nest.SetStatus(cur_gen, {"amplitude_times": t_arr[t_arr>0],
                             "amplitude_values": curr_arr[t_arr>0]})
    
    neuron = nest.Create('migliore')
    nest.SetStatus(neuron, {"V_m": -72.5})
    if stim == 'corrcost':
        nest.SetStatus(neuron,{"corrcostatratti":0,
                               "corrcost":1})
    elif stim == 'corrcostatratti':
        nest.SetStatus(neuron,{"corrcostatratti":1,
                               "corrcost":0})
    pprint.pprint(nest.GetStatus(neuron))
          
    delay = d_dt
    weight = 1.0

    conn_spec = {"rule": "all_to_all"}
    syn_spec = {'weight': weight, 'delay': delay}
    nest.Connect(cur_gen, neuron, conn_spec, syn_spec)

    multimeter = nest.Create(
        'multimeter', 1, {"interval": d_dt, 'record_from': ['V_m', 'Iadap_ini', 'Idep_ini']}) # , 'I_syn'
    nest.Connect(multimeter, neuron)
    
    tic = time.perf_counter()
    nest.Simulate(sim_time)
    toc = time.perf_counter()
    print(f"time: {toc - tic:0.4f} seconds")

    dmm = nest.GetStatus(multimeter)[0]
    V_m = dmm["events"]["V_m"]
    # I = dmm["events"]["I_syn"]
    Iada = dmm["events"]["Iadap_ini"]
    Idep = dmm["events"]["Idep_ini"]
    t = dmm["events"]["times"]

    time_stamps = np.array(t)
    membrane_voltage = np.array(V_m)
    I_adapt = np.array(Iada)
    I_dep = np.array(Idep)

    return time_stamps, membrane_voltage, I_adapt, I_dep
