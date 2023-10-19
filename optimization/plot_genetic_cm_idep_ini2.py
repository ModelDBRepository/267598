import sympy as sym
from load_eq import load_v3
import numpy as np
from geneticalgorithm import geneticalgorithm as ga
import pickle
from plotta_sim_dep_ini_Idip import plot_tr_v3_vect3_dep_ini,plot_tr_from_fit_neuron_rette_dep_ini
from load_eq import load_v3
from scipy.optimize import minimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt

def save_plot(neuron_nb,puntini,k):
    #neuron_nb='95822006'
    #k=100
    #neuron_list = []
    #puntini = False
    fitting_rule='monod'
    #puntini=False#True#puntini=True#Idep_ini

    file1 = open('dati_exp_' + neuron_nb + '.txt', 'r')
    Lines = file1.readlines()
    spk_time_orig=[]
    Istm = np.int32(Lines[3].split(','))
    #durata  e inizio stim nostre simulazioni

    durata_stim=400
    input_start_time=31.2
    limite_asse_x=431.2
    # durata  e inizio stim simulazioni Allen
    if Allen:
        durata_stim=1000
        input_start_time=1020
        limite_asse_x=1032
    tempo_finale=np.zeros(len(Istm))
    tempo_finale_tutti=np.zeros(len(Istm))
    post_tempo_finale=np.zeros(len(Istm))
    post_tempo_finale_tutti=np.zeros(len(Istm))
    corr_sp=0
    for i in range(len(Istm)):
        try:
            spk_time_orig.append(np.float64(Lines[4 + i].split(','))-input_start_time)
            corr_sp=corr_sp+1
        except:
            spk_time_orig.append([])

        if len(spk_time_orig[i])>0:
            if len(spk_time_orig[i]) == 1:
                if (durata_stim-spk_time_orig[i][len(spk_time_orig[i])-1])-(spk_time_orig[i][len(spk_time_orig[i])-1])/2>0:
                    tempo_finale[len(Istm)-i-1]=spk_time_orig[i][len(spk_time_orig[i])-1]
                    post_tempo_finale[len(Istm) - i - 1] =1.5*tempo_finale[len(Istm) - i - 1]
                tempo_finale_tutti[len(Istm) - i - 1] = spk_time_orig[i][len(spk_time_orig[i]) - 1]
                post_tempo_finale_tutti[len(Istm) - i - 1] = 1.5 * tempo_finale[len(Istm) - i - 1]
            else:
                if (durata_stim-spk_time_orig[i][len(spk_time_orig[i])-1])-(spk_time_orig[i][len(spk_time_orig[i])-1]-spk_time_orig[i][len(spk_time_orig[i])-2])*2>0:
                    tempo_finale[len(Istm)-i-1]=spk_time_orig[i][len(spk_time_orig[i])-1]
                    if len(spk_time_orig[i])>1:
                        post_tempo_finale[len(Istm)-i-1]=tempo_finale[len(Istm)-i-1]+(spk_time_orig[i][len(spk_time_orig[i])-1]-spk_time_orig[i][len(spk_time_orig[i])-2])/2
                    else:
                        post_tempo_finale[len(Istm) - i - 1] = tempo_finale[len(Istm) - i - 1] + spk_time_orig[i][len(spk_time_orig[i]) - 1] / 2

                tempo_finale_tutti[len(Istm) - i - 1] = spk_time_orig[i][len(spk_time_orig[i]) - 1]
                if len(spk_time_orig[i]) > 1:
                    post_tempo_finale_tutti[len(Istm) - i - 1] = tempo_finale_tutti[len(Istm) - i - 1] + (spk_time_orig[i][len(spk_time_orig[i]) - 1] - spk_time_orig[i][len(spk_time_orig[i]) - 2])/2
                else:
                    post_tempo_finale_tutti[len(Istm) - i - 1] = tempo_finale[len(Istm) - i - 1] + spk_time_orig[i][len(spk_time_orig[i]) - 1] / 2



        import sympy as sym
    I = sym.Symbol('I')

    I_monod_inf=-np.inf
    I_monod_sup = np.inf
    #MODIFICATO DA EMILIANO
    aux = np.logical_not(tempo_finale_tutti-tempo_finale>0)#tempo_finale > 0
    print("tempo finale",tempo_finale)
    print("tempo finale tutti", tempo_finale_tutti)

    if tempo_finale[len(Istm)-corr_sp]>0: #se la corrente meno intensa si blocca

        t_val_min= I*(post_tempo_finale_tutti[len(Istm)-corr_sp+1]-post_tempo_finale[len(Istm)-corr_sp])/(Istm[len(Istm)-corr_sp+1]-Istm[len(Istm)-corr_sp])+(Istm[len(Istm)-corr_sp+1]*post_tempo_finale[len(Istm)-corr_sp]-post_tempo_finale_tutti[len(Istm)-corr_sp+1]*Istm[len(Istm)-corr_sp])/(Istm[len(Istm)-corr_sp+1]-Istm[len(Istm)-corr_sp])

        if np.where(aux==False)[0].size>0:
            if np.where(aux==False)[0][0]>0:
                I_monod_inf=(Istm[np.where(aux==False)[0][0]-1]+Istm[np.where(aux==False)[0][0]])/2
                print("eccocci1_I_monod_inf = ", I_monod_inf)
            else:
                I_monod_inf = Istm[0] / 2
                print("eccocci2_I_monod_inf = ", I_monod_inf)

        else:
            try:
                I_monod_inf =Istm[len(Istm)-1]+(Istm[len(Istm)-1]-Istm[len(Istm)-2])/2  #se tutte le correnti sono bloccate
                a_sup=(post_tempo_finale_tutti[len(Istm)-2]-post_tempo_finale[len(Istm)-1])/(Istm[len(Istm)-2]-Istm[len(Istm)-1])
                b_sup=(Istm[len(Istm) - 2] * post_tempo_finale[len(Istm) - 1] - post_tempo_finale_tutti[len(Istm) - 2] * Istm[
                    len(Istm) - 1]) / (Istm[len(Istm) - 2] - Istm[len(Istm) - 1])
                a_inf=(post_tempo_finale_tutti[len(Istm)-corr_sp+1]-post_tempo_finale[len(Istm)-corr_sp])/(Istm[len(Istm)-corr_sp+1]-Istm[len(Istm)-corr_sp])
                b_inf=(Istm[len(Istm)-corr_sp+1]*post_tempo_finale[len(Istm)-corr_sp]-post_tempo_finale_tutti[len(Istm)-corr_sp+1]*Istm[len(Istm)-corr_sp])/(Istm[len(Istm)-corr_sp+1]-Istm[len(Istm)-corr_sp])
                print("par ", a_sup,b_sup,a_inf,b_inf)
                I_monod_inf=round((b_sup-b_inf)/(a_inf-a_sup))
                print("eccocci3_I_monod_inf = ", I_monod_inf)

                if I_monod_inf > Istm[len(Istm) - 1] or I_monod_inf < Istm[0]:
                    I_monod_inf = (Istm[np.where((tempo_finale > 0) == True)[0][len(np.where((tempo_finale > 0) == True)[0]) - 1]] + Istm[np.where((tempo_finale > 0) == True)[0][0]]) / 2
            except:
                I_monod_inf = (Istm[
                                   np.where((tempo_finale > 0) == True)[0][len(np.where((tempo_finale > 0) == True))]] +
                               Istm[np.where((tempo_finale > 0) == True)[0][0]]) / 2  # Istm[int(len(Istm) / 2)]


        #for ind in range(len(tempo_finale)-1):
        #    if tempo_finale[ind] > 0:
        #        I_monod_inf = (Istm[ind] + Istm[ind + 1]) / 2
    else:
        if tempo_finale[0]>0:
            I_monod_inf = 300

            t_val_min= I*(post_tempo_finale_tutti[1]-post_tempo_finale[0])/(Istm[1]-Istm[0])+(Istm[1]*post_tempo_finale[0]-post_tempo_finale_tutti[1]*Istm[0])/(Istm[1]-Istm[0])
            #for ind in range(len(tempo_finale)-1):
            #    if tempo_finale[ind]>0:
            #        I_monod_inf = (Istm[ind]+Istm[ind+1])/2

            if np.where(aux == False)[0].size > 0:
                if np.where(aux == False)[0][0] > 0:
                    I_monod_inf = (Istm[np.where(aux == False)[0][0]] + Istm[np.where(aux == False)[0][0] - 1]) / 2
                    print("eccocci4_I_monod_inf = ", I_monod_inf)
            else:
                I_monod_inf = Istm[len(Istm) - 1] + (Istm[len(Istm) - 1] - Istm[len(Istm) - 2]) / 2
                print("eccocci5_I_monod_inf = ", I_monod_inf)

        t_val_min = np.inf + I

    if tempo_finale[len(Istm)-1]>0:#se la corrente più intensa si blocca
        t_val_max= I*(post_tempo_finale_tutti[len(Istm)-2]-post_tempo_finale[len(Istm)-1])/(Istm[len(Istm)-2]-Istm[len(Istm)-1])+(Istm[len(Istm)-2]*post_tempo_finale[len(Istm)-1]-post_tempo_finale_tutti[len(Istm)-2]*Istm[len(Istm)-1])/(Istm[len(Istm)-2]-Istm[len(Istm)-1])#(Istm[len(Istm)-1]*post_tempo_finale[len(Istm)-2]-post_tempo_finale_tutti[len(Istm)-1]*Istm[len(Istm)-2])/(Istm[len(Istm)-2]-Istm[len(Istm)-1])
        I_monod_sup =900

        if np.where(aux==False)[0].size>0:
            if np.where(aux==False)[0][len(np.where(aux==False)[0])-1]<len(Istm)-1:#0:
                #I_monod_sup=(Istm[len(Istm)-1]+Istm[np.where(aux==False)[0][len(np.where(aux==False)[0])-1]])/2
                I_monod_sup = (Istm[np.where(aux == False)[0][len(np.where(aux == False)[0]) - 1]+1] + Istm[np.where(aux == False)[0][len(np.where(aux == False)[0]) - 1]]) / 2

        else:
            try:
                #I_monod_sup =Istm[0]/2 #se tutte le correnti sono bloccate
                a_sup = (post_tempo_finale_tutti[len(Istm) - 2] - post_tempo_finale[len(Istm) - 1]) / (
                            Istm[len(Istm) - 2] - Istm[len(Istm) - 1])
                b_sup = (Istm[len(Istm) - 2] * post_tempo_finale[len(Istm) - 1] - post_tempo_finale_tutti[
                    len(Istm) - 2] * Istm[
                             len(Istm) - 1]) / (Istm[len(Istm) - 2] - Istm[len(Istm) - 1])
                a_inf = (post_tempo_finale_tutti[len(Istm) - corr_sp + 1] - post_tempo_finale[len(Istm) - corr_sp]) / (
                            Istm[len(Istm) - corr_sp + 1] - Istm[len(Istm) - corr_sp])
                b_inf = (Istm[len(Istm) - corr_sp + 1] * post_tempo_finale[len(Istm) - corr_sp] -
                         post_tempo_finale_tutti[len(Istm) - corr_sp + 1] * Istm[len(Istm) - corr_sp]) / (
                                    Istm[len(Istm) - corr_sp + 1] - Istm[len(Istm) - corr_sp])
                I_monod_sup = round((b_sup - b_inf) / (a_inf - a_sup))
                print("eccocci3_I_monod_sup = ", I_monod_sup)
                if I_monod_sup>Istm[len(Istm)-1] or I_monod_sup<Istm[0]:
                    I_monod_sup = (Istm[np.where((tempo_finale > 0) == True)[0][len(np.where((tempo_finale > 0) == True)[0])-1]] + Istm[np.where((tempo_finale > 0) == True)[0][0]]) / 2
            except:
                I_monod_sup = (Istm[np.where((tempo_finale>0)==True)[0][len(np.where((tempo_finale>0)==True))]]+Istm[np.where((tempo_finale>0)==True)[0][0]])/2#Istm[int(len(Istm) / 2)]

        #for ind in range(len(tempo_finale)):
        #    if tempo_finale[len(tempo_finale)-ind-1] > 0:
        #        I_monod_sup = (Istm[len(tempo_finale)-ind-1] + Istm[len(tempo_finale)-ind-2]) / 2
        #        print("eccocci31_I_monod_sup = ", I_monod_sup)
    else:
        t_val_max=np.inf+I

    rette = open(neuron_nb + "_RETTE.txt", mode="w", encoding="utf-8")
    rette.write("I_monod_sup = "+str(I_monod_sup)+ "\n")
    rette.write("I_monod_inf = "+str(I_monod_inf)+ "\n")
    rette.write("t_val_max   = "+str(t_val_max)+ "\n")
    rette.write("t_val_min   = "+str(t_val_min)+ "\n")

    def exp_cum(x, a, b):
        return a * (1 - np.exp(-b * x))

    def monod(x, a, b, c, alp):
        return c + (a * np.exp(b) * x) / (alp + x)

    def exp_cum3(x, a, b, c):
        return a * (1 - np.exp(-b * x)) + c

    def monod_tot(x):
        err = []
        for i in range(corr_sp2):
            err.append(sum((x[2] + ((x[0] * np.exp(x[1] * Istm[len(Istm) - 1 - i] / 1000).astype(np.float64) * np.array(t_sp_abs_tutti[corr_sp2 - i - 1]).astype(np.float64))) / (x[3] + np.array(t_sp_abs_tutti[corr_sp2 - i - 1]).astype(np.float64)) - np.array(Iada_tutti[corr_sp2 - i - 1]).astype(np.float64)) ** 2))  # *np.array(t_var_tutti[corr_sp2-i-1]).astype(np.float64)))
        err_tot = sum(err)
        return err_tot.astype(np.float64)


    with open('neuron_' + neuron_nb + '_cm_idep.pkl', 'rb') as f:  # Python 3: open(..., 'wb')

        [t_sp_abs_tutti, Iada_tutti, Istm_c, spk_time_c, t_sp_sim_c, Vconvfact, vth, vrm, bet, delta1, sc,time_scale, Idep_ini_vr, input_start_time, dur_sign, fitta, popt,corr_sp_c,corr_sp2_c,tao_m,Ith,time_soglia,Idep_ini] = pickle.load(f)
        #Idep_ini=0.9
        
    if puntini:
        t_sp_sim_c = []
        #ada_sim_c = []
        tspkIadap = open(neuron_nb + "_TSPK_IADAP.txt", mode="w", encoding="utf-8")
        tspkIadap_final=open(neuron_nb + "_TSPK_IADAP_FINAL.txt", mode="w", encoding="utf-8")
        for i in range(corr_sp_c):  #
            aux = Istm_c / sc
            alpha = aux.tolist()
            if i<corr_sp_c-corr_sp2_c:
                [t_aux, time, voltage,Iada_final] = plot_tr_v3_vect3_dep_ini(Vconvfact, vth, vrm, [], bet, delta1,alpha[len(alpha) - corr_sp_c + i], sc, time_scale, Idep_ini_vr,input_start_time, dur_sign,Idep_ini,Ith)
            else:
                [t_aux, time, voltage,Iada_final] = plot_tr_v3_vect3_dep_ini(Vconvfact, vth, vrm, Iada_tutti[i-(corr_sp_c-corr_sp2_c)], bet, delta1,alpha[len(alpha) - corr_sp_c + i], sc, time_scale, Idep_ini_vr,input_start_time, dur_sign,Idep_ini,Ith)
                Iaaux = Iada_tutti[i - (corr_sp_c - corr_sp2_c)]
                for j in range(len(Iaaux)):

                    try:
                        tspkIadap.write(str(Istm_c[len(alpha) - corr_sp_c + i]) + "  " + str(t_aux[j]+2) + "  " + str(Iaaux[j]) + "\n")
                    except:
                        print(neuron_nb)
            
            for j in range(len(Iada_final)):
                tspkIadap_final.write(str(Istm_c[len(alpha) - corr_sp_c + i]) + "  " + str(t_aux[j]) + "  " + str(Iada_final[j]) + "\n")

                
            t_sp_sim_c.append(t_aux)
            #ada_sim_c.append(ada_aux)
    else:
        t_sp_sim_c=[]
        ada_sim_c=[]
        xdata = np.linspace(0, int(input_start_time + dur_sign) + 30, 100)
        with open('best_res_' + neuron_nb + '_cm.pkl', 'rb') as f:
            [best_res] = pickle.load(f)

        tspkmonod = open("MONOD_"+neuron_nb + "_TSPK.txt", mode="w", encoding="utf-8")
        for i in range(corr_sp_c):  #

            if fitting_rule=='monod':
                a = [0, 0, -200, 0]
                b = [20000, 10, 200, 100000]
                fun_loss_sel = monod_tot
                fun_sel = monod
                popt[0] = best_res.x[0]
                popt[1] = best_res.x[1] * Istm_c[len(Istm_c)-corr_sp_c+i]/1000
                popt[2] = best_res.x[2]
                popt[3] = best_res.x[3]

                aux = Istm_c / sc
                alpha=aux.tolist()
                is_Neuron_tr = True

                [t_aux,ada_aux,time,voltage]=plot_tr_from_fit_neuron_rette_dep_ini(Vconvfact, vth, vrm, fun_sel , popt, bet, delta1, alpha[len(alpha)-corr_sp_c+i],sc, time_scale, Idep_ini_vr, input_start_time, dur_sign,t_val_min,I_monod_inf,t_val_max,I_monod_sup,Idep_ini,Ith)#np.inf
                t_sp_sim_c.append(t_aux)
                ada_sim_c.append(ada_aux)
                #mi salvo i tempi degli spikes della MONOD
                for j in range(len(t_aux)):
                    tspkmonod.write(str(Istm_c[len(alpha) - corr_sp_c + i]) + "  " + str(t_aux[j]) + "\n")

                plt.figure(k+3);
                if popt.size == 2:
                    plt.plot(xdata, fun_sel(xdata, *popt), label='fit: a=%5.3f, b=%5.3f' % tuple(popt))

                if popt.size == 3:
                    plt.plot(xdata, fun_sel(xdata, *popt), label='fit: a=%5.3f, b=%5.3f,c=%5.3f' % tuple(popt))

                if popt.size == 4:
                    plt.plot(xdata, fun_sel(xdata, *popt), label='fit: a=%5.3f, b=%5.3f,c=%5.3f,alpha=%5.3f' % tuple(popt))

                plt.legend(bbox_to_anchor=(0.2, 0.2), loc='upper left', borderaxespad=0.)
                plt.xlabel('Time (ms)')
                plt.ylabel('Iadap')
                
    tspkpuntini = open("PUNTINI_"+neuron_nb + "_TSPK.txt", mode="w", encoding="utf-8")
    for i in range(corr_sp_c):
    
        plt.figure(k+1);
        #plt.scatter(spk_time_c[i], Istm_c[len(Istm_c)-1-i] * np.ones(len(spk_time_c[i])), marker='|',color='r');
        plt.scatter(spk_time_orig[i], Istm_c[len(Istm_c) - 1 - i] * np.ones(len(spk_time_orig[i])), marker='|', color='r');
        plt.scatter(t_sp_sim_c[len(t_sp_sim_c)-1-i],Istm_c[len(Istm_c)-1-i] * np.ones(len(t_sp_sim_c[len(t_sp_sim_c)-1-i])), marker='|', color='b');
        plt.xlim(0,limite_asse_x)
        plt.title('neuron ' + neuron_nb)
        plt.ylabel('Current (pA)')
        plt.xlabel('spike times (ms)')
        if puntini:
            plt.savefig('raster_plot_' + neuron_nb + '.png')
            #plt.savefig('sper_raster_plot_' + neuron_nb + '.png')
        else:
            plt.savefig('MONOD_raster_plot_' + neuron_nb + '.png')

        plt.figure(k+2)
        plt.scatter(Istm_c[len(Istm_c)-1-i],len(spk_time_c[i]),marker='*',color='r')
        plt.scatter(Istm_c[len(Istm_c)-1-i],len(t_sp_sim_c[len(t_sp_sim_c)-1-i]),marker='*',color='b')
        plt.title('neuron ' + neuron_nb)
        plt.ylabel('number of spikes')
        plt.xlabel('Current (pA)')
        if puntini:
            plt.savefig('spike_numbers_' + neuron_nb + '.png')
            #plt.savefig('sper_spike_numbers_' + neuron_nb + '.png')
            #mi salvo i tempi degli spikes dei PUNTINI
            for j in range(len(t_sp_sim_c[i])):
                tspkpuntini.write(str(Istm_c[len(alpha) - corr_sp_c + i]) + "  " + str(t_sp_sim_c[i][j]) + "\n")
        else:
            plt.savefig('MONOD_spike_numbers_' + neuron_nb + '.png')

    for i in range(corr_sp2_c):
        plt.figure(k+3);
        plt.scatter(t_sp_abs_tutti[i], Iada_tutti[i])
        plt.title('neuron ' + neuron_nb)
        plt.xlabel('Time (ms)')
        plt.ylabel('Iadap')
        plt.xlim(0,limite_asse_x)
        #plt.title('Iada_tot Idep' + str(Idep_ini_vr) + ' tao ' + str(tao_m) + ' Ith ' + str(Ith))
        if puntini:
            plt.savefig('punti_' + neuron_nb + '.png')
        else:
            plt.savefig('MONOD_punti_' + neuron_nb + '.png')
            
neuron_list=['96711008','99111002','99111001','99111000','97911000','97911001','97911002','97428000','97428001','98205021','98205022','98205024','98205025','97509008','97509009','97509010','97509011','97717005','98D15008','98D15009','99111004','99111006','98513011','95817000','95817001','95817002','95810006','95810006','95810007','95810008','95810010','95810011','95810022','95810023','95810024','95810025','95810030','95810031','95810032','95810033','95810037','95810038','95810039','95810040','95810041','95817003','95817004','95817005','95817006','95817007','95817008','95822002','95822003','95822005','95822006','95822010','95822011','95824000','95824004','95824006','95831000','95831001','95831002','95831003','95831004','95912004','95912005','95912006','95912007','95914001','95914002','95914003','95914004','95810012','95810013','95810014','95810015','95810026','95810027','95810028','95810029','95822000','95822001','95822009']

Allen=False
if Allen:
    neuron_list=['47','48']

puntini=False
k=1
for i in neuron_list:
    save_plot(i, puntini,k)
    k=k+10
