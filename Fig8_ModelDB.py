import numpy as np
import matplotlib.pyplot as plt
from mds_python_model import mds_sim

fig = plt.figure()
p_time_stamps = []
p_membrane_voltage = []
n_time_stamps = []
n_membrane_voltage = []
Istim = ['A','B','C','D','E','F']
sim_length=1000
campionamento=40
d_dt=0.005*campionamento
for i in range(1,7):
    if Istim[i-1]=='A':
        Istim0=0
        Istim1=600
        Istim2=400
        Istim3=1000
        current = np.ones(int(sim_length/d_dt))*Istim0
        change_cur = 200
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim1
        change_cur = 300
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim2
        change_cur = 500
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim3
        change_cur = 600
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim0
        change_cur = 800
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim2
    if Istim[i-1]=='A':
        Istim0=0
        Istim1=600
        Istim2=400
        Istim3=1000
        current = np.ones(int(sim_length/d_dt))*Istim0
        change_cur = 200
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim1
        change_cur = 300
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim2
        change_cur = 500
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim3
        change_cur = 600
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim0
        change_cur = 800
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim2
    if Istim[i-1]=='B':
        Istim0=0
        Istim1=400
        Istim2=700
        Istim3=200
        Istim4=1000
        current = np.ones(int(sim_length/d_dt))*Istim0
        change_cur = 200
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim1
        change_cur = 400
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim2
        change_cur = 600
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim3
        change_cur = 800
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim4
    if Istim[i-1]=='C':
        Istim0=0
        Istim1=600
        Istim2=500
        Istim3=250
        Istim4=1000
        current = np.ones(int(sim_length/d_dt))*Istim0
        change_cur = 200
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim1
        change_cur = 400
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim2
        change_cur = 600
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim3
        change_cur = 800
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim4
    if Istim[i-1]=='D':
        Istim0=0
        Istim1=600
        Istim2=500
        Istim3=800
        Istim4=1000
        current = np.ones(int(sim_length/d_dt))*Istim0
        change_cur = 200
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim1
        change_cur = 400
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim2
        change_cur = 600
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim3
        change_cur = 800
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim4
    if Istim[i-1]=='E':
        Istim0=0
        Istim1=600
        Istim2=500
        Istim3=400
        Istim4=800
        Istim5=1000
        current = np.ones(int(sim_length/d_dt))*Istim0
        change_cur = 200
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim1
        change_cur = 240
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim2
        change_cur = 300
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim1
        change_cur = 340
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim3
        change_cur = 400
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim4
        change_cur = 800
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim5
    if Istim[i-1]=='F':
        Istim0=0
        Istim1=1000
        Istim2=800
        Istim3=400
        current = np.ones(int(sim_length/d_dt))*Istim0
        change_cur = 200
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim1
        change_cur = 400
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim2
        change_cur = 600
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim3
        change_cur = 800
        current[int(change_cur/d_dt):int(sim_length/d_dt) +1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim1
    plt.subplot(6,1,i)
    a, b , p_I_adapt, p_I_dep, p_monod_plot, p_Iadap0max_plot, p_init_sign_plot = mds_sim('corrcostatratti',current, str(Istim[i-1]),d_dt)
    p_time_stamps.append(a)
    p_membrane_voltage.append(b)

    plt.subplot(6,1,i)
    plt.ylabel(str(Istim[i-1]),fontsize=15,rotation='horizontal', ha='right',va="center",weight='bold')
    plt.plot(p_time_stamps[i-1], p_membrane_voltage[i-1], 'k', label='python')
    plt.xlim([0, 1000])
    if i<6:
        plt.xticks(color='white')
    plt.ylim([-75, -30])
    plt.xlabel('Time (ms)')
#plt.show()
plt.savefig('Model_traces_for_piecewise_currents.png')
