import numpy as np
import matplotlib.pyplot as plt
from mds_python_model import mds_sim

fig = plt.figure()
p_time_stamps = []
p_membrane_voltage = []
n_time_stamps = []
n_membrane_voltage = []
Istim = [400, 500, 600, 700, 800, 1000, 1100]
sim_length=1000
campionamento=40
d_dt=0.005*campionamento
for i in range(1,8):
    Istim0 = 0
    current = np.ones(int(sim_length/d_dt))*Istim0
    change_cur = 400
    current[int(change_cur/d_dt):int(sim_length/d_dt) +
        1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim[i-1]
    change_cur = 800
    current[int(change_cur/d_dt):int(sim_length/d_dt) +
        1] = np.ones(len(current[int(change_cur/d_dt):int(sim_length/d_dt)+1]))*Istim0
    plt.subplot(7,1,i)
    a, b , p_I_adapt, p_I_dep, p_monod_plot, p_Iadap0max_plot, p_init_sign_plot = mds_sim('corrcost',current, str(Istim[i-1]),d_dt)
    p_time_stamps.append(a)
    p_membrane_voltage.append(b)

    plt.subplot(7,1,i)
    plt.xlim([0, 1000])
    if i<7:
        plt.xticks(color='white')
    plt.plot(p_time_stamps[i-1], p_membrane_voltage[i-1], 'k', label='python')
    plt.xlabel('Time (ms)')
    plt.ylabel(str(Istim[i-1]/1000)+' nA',fontsize=7,rotation='horizontal', ha='right',va="center",weight='bold')
plt.ylim([-75, -30])
#plt.show()
plt.savefig('Model_traces_for_constant_current_injections.png')
