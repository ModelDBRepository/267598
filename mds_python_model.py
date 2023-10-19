import matplotlib.pyplot as plt
import numpy as np
import time
import sys

def tagliorette(corr,sim_type):
    dur_sign = np.inf
    vinc_inf = 700

    if corr < vinc_inf and corr > 0:
        dur_sign = lin_func_inf = 0.68*corr - 190.0
    if sim_type=='corrcostatratti':
        vinc_sup = np.inf

        if corr > vinc_sup:
            dur_sign = lin_func_sup = corr + np.inf
    elif sim_type == 'corrcost':
        vinc_sup = 1300

        if corr > vinc_sup:
            dur_sign = lin_func_sup = 76.86-0.028*corr
    else:
        print("#####################!!! Sim_type not defined!!!!#######")
        return
    
    return dur_sign


def V(t, delta, Psi, alpha, beta, IaA0, IdA0, t0, V0):
    t_step = t - t0
    VV_1 = 0.5 / ((beta -delta) * (pow(beta,2) + (beta-1.0) * delta) * (4.0 * beta - (1.0 + delta) ** 2.0)) * Psi
    VV_2 = 2.0 * np.exp(-t_step * beta) * (beta-1.0) * beta * (beta - delta) * Psi
    VV_3 = (pow(beta,2) + ((- 1.0) + beta) * delta) * Psi
    VV_4 = np.exp((1.0 / 2.0) * t_step * (-1.0 + delta -Psi))
    VV_5 = beta * (beta -delta) * (-1.0 -delta + beta * (3.0 + delta -Psi) + Psi)
    VV_6 = (pow(beta,2) -delta + beta * delta)
    VV_7 = (1.0 + (-2.0) * beta + delta -Psi)
    VV_8 = np.exp((1.0 / 2.0) * t_step * (-1.0 + delta + Psi))
    VV_9 = beta * (beta-delta) * (-1.0 -delta -Psi + beta * (3.0 + delta + Psi))
    VV_10 = (pow(beta,2) - delta + beta * delta)
    VV_11 = (1.0 + (-2.0) * beta + delta + Psi)
    return VV_1 * (VV_2 * IdA0 + -2.0 * (alpha -beta + delta) * VV_3 + VV_4 * (IdA0 * VV_5 - VV_6 * (alpha * VV_7 + (beta -delta) * (-1.0 + 2.0 * IaA0 * beta -delta + Psi + V0 * (-1.0 -delta + Psi)))) + VV_8 * (-IdA0 * VV_9 + VV_10 * (alpha * VV_11 + (beta -delta) * (-1.0 + 2.0 * IaA0 * beta-delta -Psi -V0 * (1.0 + delta + Psi)))))


def Iadap(t, delta, Psi, alpha, beta, IaA0, IdA0, t0, V0):
    AA_1 = -4.0*beta**3.0+beta**2.0*(-1.0+delta)**2.0-delta*(1.0+delta)**2.0+beta*delta*(5+2.0*delta+delta**2.0)
    AA_2 = 2.0*np.exp(-(t - t0) * beta)*beta*(4.0*beta**2.0+delta*(1.0+delta)**2.0-beta*(1.0+6*delta+delta**2.0))
    AA_3 = np.exp((1.0 / 2.0)*(t-t0)*(-1.0+delta+Psi))
    AA_4 = -1.0-2.0*delta-delta**2.0-Psi+delta*Psi+2.0*beta*(2.0+Psi)
    AA_5 = (beta**2.0-delta+beta*delta)
    AA_6 = (1.0+(-4.0)*beta+2.0*delta+delta**2.0+Psi-delta*Psi)
    AA_7 = (1.0+delta)*(-1.0-delta+Psi)
    AA_8 = np.exp((-1.0)*(1.0 / 2.0) * (t-t0) * (1.0-delta+Psi))
    AA_9 = beta * (beta-delta)*(1.0+2.0*delta+delta**2.0-Psi+delta*Psi+2.0*beta*(-2.0+Psi))
    AA_10 = (beta**2.0-delta+beta*delta)
    AA_11 = (1.0-4.0*beta+2.0*delta+delta**2.0-Psi+delta*Psi)
    AA_12 = (2.0*(beta-delta)*(beta**2.0+(-1.0+beta)*delta)*(4.0*beta-(1.0+delta)**2.0))
    to_return = (-2.0 * alpha * AA_1 + IdA0 * AA_2 + AA_3 * (-IdA0 * beta * (beta - delta) * AA_4 + AA_5 * (alpha * AA_6 + (beta - delta) * (4.0 * IaA0 * beta -2.0 * (1.0 + V0) * Psi + IaA0 * AA_7))) + AA_8 * (IdA0 * AA_9 + AA_10 * (alpha * AA_11 - (beta-delta) * (-4.0 * IaA0 * beta - 2.0*(1.0+V0)*Psi+IaA0*(1.0+delta)*(1.0+delta+Psi)))))/ AA_12
    return to_return

def Idep(t, beta, IdA0, t0):
    return np.exp(((-1) * t + t0) * beta) * IdA0


def exp_cum(x, a, b):
    return a * (1 - np.exp(-b * x))


def monod(x, a, b, c, alp):
    to_return = c + (a * np.exp(b) * x) / (alp + x)
    return to_return

def mds_sim(sim_type, total_current, sim_id, d_dt):
    EL = -72.5
    vres = -65
    vtm = -52
    Cm = 2047.4164432004916
    ith = 300.9987901902274
    tao_m = 3310.462136574088
    sc = 4526.328798037026
    bet = 0.24522251335200956
    delta1 = 0.009906254244852036
    cost_idep_ini = 0.017625482908326662
    Idep_ini_vr = 1.0122905259090516
    psi1 = 0.1975362978159442
    a=14.2787
    b=-2.10966
    c=0.0608809
    alp=225.491
    istim_min_spikinig_exp=400
    istim_max_spikinig_exp=1000
    time_scale = 1 / (-sc / (Cm * EL))
    H = (90+EL)*sc*(bet-delta1)/(EL*(-200))
    cor = total_current
    sim_length = 1000
    if sim_type=='corrcostatratti':
        corrcostatratti = 1
        corrcost = 0
    elif sim_type == 'corrcost':
        corrcostatratti = 0
        corrcost = 1
    else:
        print("#####################!!! Sim_type not defined!!!!#######")
        return
    tic = time.perf_counter()

    Vconvfact = -EL
    vth = vtm/Vconvfact
    vrm = vres/Vconvfact

    t0_val = 0
    vini_neg = EL

    ts = np.inf

    dt = d_dt/time_scale
    init_sign = 0
    ref_t = 2

    t0_val = 0
    psi1 = ((-4)*bet+((1+delta1)**2))**(0.5)

    Idep_ini = 0
    Iadap_ini = 0
    out = []
    t_out = []

    t_final = t0_val+dt
    v_ini = -1
    vini_prec = v_ini
    mul = 15

    f = open('t_spk_simulated_SIM4before_contheta_newinit_codicepulito'+sim_id+'.txt', 'w')
    i = 0

    soglia_sign = 10
    Ide = []
    Iada = []
    Ide2 = []
    Iada2 = []
    tetalist = []
    monod_plot = []
    Iadap0max_plot = []
    init_sign_plot = []
    
    t_spk = -3*d_dt
    afirst = 0
    meancorlastis = 0
    stdcorlastis = 0
    sis = 0

    firstSpikeFlag = False
    counter = 0

    while(t_final*time_scale < sim_length):
        if (t_final-init_sign)*time_scale >= tagliorette(cor[i],sim_type):
            if corrcostatratti:
                if cor[i] > ith:
                    if cor[i-1] < ith or i == 0:
                        init_sign = t_final
                        Idep_ini = cost_idep_ini*(cor[i]-ith)
                        Iadap_ini = 0
                    if cor[i-1] > ith and cor[i-1] < cor[i]:
                        init_sign = init_sign*(1+(cor[i-1]-ith)/cor[i-1])
                        Idep_ini = cost_idep_ini*(cor[i]-ith)
            if corrcost:
                if cor[i] > ith:
                    if cor[i-1] < ith or i == 0:
                        init_sign = t_final
                        Idep_ini = cost_idep_ini*(cor[i]-ith)
                        Iadap_ini = 0
            if cor[i-1] == 0:
                v_ini = vini_prec
            else:
                v_ini= (EL+(1-np.exp(-cor[i]/1000))*(vtm-EL))/Vconvfact
            vini_prec = v_ini
        else:
            vini_prec = v_ini
            if (cor[i] < ith and cor[i] >= 0) or i == 0:
                if cor[i-1] < 0:
                    Iadap_ini = 90/EL + 1
                    Idep_ini = 0
                    v_ini = vini_prec
                if ((cor[i] / sc) / (bet - delta1) - 1) <= v_ini:
                    Idep_ini = 0
                    Iadap_ini = (cor[i] / sc) / (bet - delta1)
                    v_ini = ((cor[i] / sc) / (bet - delta1) - 1)
                else:
                    v_ini = V(t_final, delta1, psi1,
                              cor[i]/sc, bet, Iadap_ini, Idep_ini, t0_val, v_ini)
                    Iadap_ini = Iadap(
                        t_final, delta1, psi1, cor[i] / sc, bet, Iadap_ini, Idep_ini, t0_val, v_ini)
                    Idep_ini = Idep(t_final, bet, Idep_ini, t0_val)
                if v_ini * Vconvfact < -90:
                    v_ini = -90 / Vconvfact
                    Iadap_ini = 0
            else:
                if cor[i] < cor[i-1] and cor[i] > 0 and (t_spk+2*d_dt) < t_final*time_scale:
                    teta = (out[i-1]/(cor[i-1] / sc))*(1/dt-delta1) - \
                        (out[i-2]/((cor[i-1] / sc)*dt))-delta1/(cor[i-1] / sc)-1
                    if teta < 0:
                        teta = 0
                    Idep_ini = Iadap_ini + teta * (cor[i] / sc) / bet
                    tetalist.append(teta)
                    v_ini = V(t=t_final, delta=delta1, Psi=psi1,
                              alpha=cor[i]/sc, beta=bet, IaA0=Iadap_ini, IdA0=Idep_ini, t0=t0_val, V0=v_ini)
                    Iadap_ini = Iadap(t=t_final, delta=delta1, Psi=psi1,
                                      alpha=cor[i]/sc, beta=bet, IaA0=Iadap_ini, IdA0=Idep_ini, t0=t0_val, V0=v_ini)
                    Idep_ini = Idep(t=t_final, beta=bet,
                                    IdA0=Idep_ini, t0=t0_val)
                else:
                    if cor[i] > 0:
                        v_ini = V(t=t_final, delta=delta1, Psi=psi1,
                                  alpha=cor[i]/sc, beta=bet, IaA0=Iadap_ini, IdA0=Idep_ini, t0=t0_val, V0=v_ini)
                        Iadap_ini = Iadap(t=t_final, delta=delta1, Psi=psi1,
                                          alpha=cor[i]/sc, beta=bet, IaA0=Iadap_ini, IdA0=Idep_ini, t0=t0_val, V0=v_ini)
                        Idep_ini = Idep(t=t_final, beta=bet,
                                        IdA0=Idep_ini, t0=t0_val)
                if cor[i-1] != cor[i] and (cor[i] < 0 and cor[i] > -200):
                    Iadap_ini = (90+EL)*cor[i]/(EL*(-200))
                    Idep_ini = 0
                    v_ini = vini_prec
                if cor[i] < 0 and cor[i] > -200:
                    v_ini = V(t=t_final, delta=delta1, Psi=psi1, alpha=H *
                              cor[i]/sc, beta=bet, IaA0=Iadap_ini, IdA0=Idep_ini, t0=t0_val, V0=v_ini)
                    Iadap_ini = Iadap(t=t_final, delta=delta1, Psi=psi1, alpha=H *
                                      cor[i]/sc, beta=bet, IaA0=Iadap_ini, IdA0=Idep_ini, t0=t0_val, V0=v_ini)
                    Idep_ini = Idep(t=t_final, beta=bet,
                                    IdA0=Idep_ini, t0=t0_val)
                if cor[i-1] != cor[i] and cor[i] <= -200:
                    Iadap_ini = 90/EL + 1
                    Idep_ini = 0
                    v_ini = vini_prec
                if cor[i] <= -200:
                    v_ini = V(t=t_final, delta=delta1, Psi=psi1, alpha=H *
                              cor[i]/sc, beta=bet, IaA0=Iadap_ini, IdA0=Idep_ini, t0=t0_val, V0=v_ini)
                    Iadap_ini = Iadap(t=t_final, delta=delta1, Psi=psi1, alpha=H *
                                      cor[i]/sc, beta=bet, IaA0=Iadap_ini, IdA0=Idep_ini, t0=t0_val, V0=v_ini)
                    Idep_ini = Idep(t=t_final, beta=bet,
                                    IdA0=Idep_ini, t0=t0_val)
                if v_ini*Vconvfact < -90:
                    v_ini = -90/Vconvfact
                    Iadap_ini = 0

            if corrcostatratti:
                if cor[i] > ith:
                    if cor[i-1] < ith or i == 0:
                        init_sign = t_final
                        Idep_ini = cost_idep_ini*(cor[i]-ith)
                        Iadap_ini = 0
                    if cor[i-1] > ith and cor[i-1] < cor[i]:
                        init_sign = init_sign*(1+(cor[i-1]-ith)/cor[i-1])
                        Idep_ini = cost_idep_ini*(cor[i]-ith)
                        Iadap_ini = 0
            if corrcost:
                if cor[i] > ith:
                    if cor[i-1] < ith or i == 0:
                        init_sign = t_final
                        Idep_ini = cost_idep_ini*(cor[i]-ith)
                        Iadap_ini = 0
            
            if v_ini > vth:
                t_spk = t_final*time_scale
                f.write(str(round(t_spk, 3)) + ' \n')
                v_ini = vrm

                print('***spike***')
                print('t ', t_final, 'val_ist V', v_ini * Vconvfact, 'adap',
                      Iadap_ini, 'adap', Iadap_ini, 't_ini', init_sign)
                print('************')
                if cor[i] < istim_min_spikinig_exp or cor[i] > istim_max_spikinig_exp:
                    if corrcost or corrcostatratti:
                        if firstSpikeFlag == False or cor[i-1]!=cor[i]:
                            firstSpikeFlag=True
                            paramL = monod((t_final-init_sign)*time_scale,a,b*cor[i]/1000,c,alp)
                            if paramL<0:
                                if a > 0:
                                    c_aux = c - paramL
                                else:
                                    c_aux = -a*np.exp(b*cor[i]/1000)
                            else:
                                c_aux = c
                        Iadap_ini = monod((t_final-init_sign) *
                              time_scale, a, b*cor[i]/1000, c_aux, alp)

                    else:
                        c_aux = c
                        Iadap_ini = monod((t_final-init_sign) *
                              time_scale, a, b*cor[i]/1000, c_aux, alp)
                else:
                    Iadap_ini = monod((t_final-init_sign) *
                                      time_scale, a, b*cor[i]/1000, c, alp)
                    if Iadap_ini<0:
                        if firstSpikeFlag == False or cor[i-1]!=cor[i]:
                            firstSpikeFlag=True
                            paramL = monod((t_final-init_sign)*time_scale,a,b*cor[i]/1000,c,alp)                      
                        if counter<2:
                            Iadap_ini = 0
                            counter = counter + 1
                        elif counter == 2:
                            c_aux = c
                            if a > 0:
                                c_aux = c - paramL
                            else:
                                c_aux = -a*np.exp(b*cor[i]/1000)
                            Iadap_ini = monod((t_final-init_sign) * time_scale, a, b*cor[i]/1000, c_aux, alp)
                            counter = counter + 1
                        else:
                            Iadap_ini = monod((t_final-init_sign) * time_scale, a, b*cor[i]/1000, c_aux, alp)
                if cor[i] < ith:
                    Idep_ini = 0
                    Iadap_ini = 0
                else:
                    Idep_ini = Idep_ini_vr
                for k in range(int(ref_t / d_dt)):
                    out.append(v_ini)
                    t_out.append(t_final)
                    t_final = t_final + dt
                    Iada.append(Iadap_ini)
                    Ide.append(Idep_ini)
                    
                    i = i + 1
            vini_prec = v_ini

        out.append(v_ini)
        t_out.append(t_final)
        Iada.append(Iadap_ini)
        Ide.append(Idep_ini)
        i = i + 1
        t0_val = t_final
        t_final = t0_val+dt

    membrane_voltage = np.array(out) * Vconvfact
    I_adapt = np.array(Iada)
    I_dep = np.array(Ide)
    monod_plot = np.array(monod_plot)
    Iadap0max_plot = np.array(Iadap0max_plot)
    init_sign_plot = np.array(init_sign_plot)* time_scale

    return np.array(t_out)*time_scale, membrane_voltage, I_adapt, I_dep, monod_plot, Iadap0max_plot, init_sign_plot
