import sympy as sym
from load_eq import load_v3
import numpy as np
from geneticalgorithm import geneticalgorithm as ga
import pickle
from plotta_sim_dep_ini import plot_tr_v3_vect3_dep_ini2,plot_tr_from_fit_neuron_dep_ini2
from load_eq import load_v3
from scipy.optimize import minimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt

def runsimulation(neuron_nb,cm_fixed):

    def exp_cum(x, a, b):
        return a * (1 - np.exp(-b * x))


    def monod(x, a, b, c, alp):
        return c + (a * np.exp(b) * x) / (alp + x)


    def exp_cum3(x, a, b, c):
        return a * (1 - np.exp(-b * x)) + c


    def pol2(x, a, b):
        return a * x + b * x ** 2


    def monod_tot(x):
        err = []
        for i in range(corr_sp2):
            err.append(sum((x[2] + ((x[0] * np.exp(x[1] * Istm[len(Istm) - 1 - i] / 1000).astype(np.float64) * np.array(t_sp_abs_tutti[corr_sp2 - i - 1]).astype(np.float64))) / (x[3] + np.array(t_sp_abs_tutti[corr_sp2 - i - 1]).astype(np.float64)) - np.array(Iada_tutti[corr_sp2 - i - 1]).astype(np.float64)) ** 2))  # *np.array(t_var_tutti[corr_sp2-i-1]).astype(np.float64)))
        err_tot = sum(err)
        return err_tot.astype(np.float64)



    def ida_to_mat(IdA_ex):
        IdA_ex_tras=(IdA_ex+(-IdA_min))/(IdA_max-IdA_min)*(n-1)
        return IdA_ex_tras

    def iaa_to_mat(IaA_ex):
        IaA_ex_tras=(IaA_ex+(-IaA_min))/(IaA_max-IaA_min)*(m-1)
        return IaA_ex_tras

    def mat_to_ida(IdA_ex_tras):
        IdA_ex=IdA_min+IdA_ex_tras*(IdA_max-IdA_min)/(n-1)
        return IdA_ex

    def mat_to_iaa(IaA_ex_tras):
        IaA_ex=IaA_min+IaA_ex_tras*(IaA_max-IaA_min)/(m-1)
        return IaA_ex


    def loss_func(x):
        par_sc=x[0]
        tao_m=x[1]
        Ith=x[2]
        Idep_ini_vr=x[3]
        Cm =x[4] #189.79
        Idep_ini=x[5]

        t = sym.Symbol('t')
        delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')
        [V, Iadap, Idep] = load_v3()
        min_sc = (-Cm * EL) / tao_m + (2 / (1 + vth)) * (Ith + np.sqrt(Ith * (Ith - Cm * EL * (1 + vth) / tao_m)))
        sc = min_sc + par_sc
        k_adap_min = sc * (-EL / tao_m + Ith / (Cm * (1 + vth))) / (EL ** 2)
        k_adap_max = ((Cm * EL - sc * tao_m) ** 2) / (4 * Cm * (EL ** 2) * (tao_m ** 2))
        k_adap = k_adap_min + (k_adap_max - k_adap_min) * 0.01  #
        delta1 = -Cm * EL / (sc * tao_m)
        bet = Cm * (EL ** 2) * k_adap / (sc ** 2)
        psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)
        time_scale = 1 / (-sc / (Cm * EL))
        Vconvfact=-EL
        aux = Istm / sc
        alpha3 = aux.tolist()
        spt_00 = np.zeros(len(alpha3))
        spt_max = np.zeros(len(alpha3))
        spt_min = np.zeros(len(alpha3))
        err_fs= np.zeros(len(alpha3))
        err_ss_min= np.zeros(len(alpha3))
        err_ss_max= np.zeros(len(alpha3))

        for i in range(corr_sp):
            aux = V.subs(alpha, alpha3[len(alpha3) - 1 - i]).subs(beta, bet).subs(delta, delta1).subs(t0, 0).subs(V0,-1).subs(IaA0, 0).subs(IdA0, Idep_ini*(Istm[len(alpha3) - 1 - i]-Ith)*(Istm[len(alpha3) - 1 - i]>Ith)).subs(Psi, psi1)
            lam_x = sym.lambdify(t, aux, modules=['numpy'])
            x_vals = np.linspace(0, 1000, 10001) / time_scale
            y_vals = lam_x(x_vals)
            aus = np.nonzero((y_vals  > vth) * y_vals)
            if aus[0].size > 0:
                spt_00[i] = x_vals[aus[0][0]] * time_scale
            else:
                spt_00[i] = -10000

        for i in range(corr_sp2):
            IaA_max = Idep_ini_vr+ alpha3[len(alpha3) - 1 - i] / bet + (delta1 / bet) * (1 + vrm)
            aux = V.subs(alpha, alpha3[len(alpha3) - 1 - i]).subs(beta, bet).subs(delta, delta1).subs(t0, 0).subs(V0, vrm).subs(IaA0, IaA_max).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)
            lam_x = sym.lambdify(t, aux, modules=['numpy'])
            x_vals = np.linspace(0, int(2*max(spk_interval[i][1:])), int(2*(max(spk_interval[i][1:])))+1) / time_scale
            y_vals = lam_x(x_vals)
            aus = np.nonzero((y_vals  > vth) * y_vals)
            if aus[0].size > 0:
                spt_max[i] = x_vals[aus[0][0]] * time_scale
            else:
                spt_max[i] = 10000


        for i in range(corr_sp2):
            IaA_min=0
            aux = V.subs(alpha, alpha3[len(alpha3) - 1 - i]).subs(beta, bet).subs(delta, delta1).subs(t0, 0).subs(V0,vrm).subs(IaA0, IaA_min).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)
            lam_x = sym.lambdify(t, aux, modules=['numpy'])
            x_vals = np.linspace(0,int(2*min(spk_interval[i][1:])),int(2*min(spk_interval[i][1:]))+1 ) / time_scale
            y_vals = lam_x(x_vals)
            aus = np.nonzero((y_vals > vth) * y_vals)
            if aus[0].size > 0:
                spt_min[i] = x_vals[aus[0][0]] * time_scale
            else:
                spt_min[i] = 10000

            for i in range(corr_sp):

                err_fs[i]=abs(spt_00[i]-spk_interval[i][0])
                if (i<corr_sp2):
                    err_ss_min[i] = sum(np.maximum(0, spt_min[i]-spk_interval[i][1:]))#max(0,(spt_min[i] - min(spk_interval[i][1:])))
                    err_ss_max[i] = sum(np.maximum(0, spk_interval[i][1:]-spt_max[i]))#max(0,max(spk_interval[i][1:]-spt_max[i]))


                else:
                    err_ss_min[i] = 0
                    err_ss_max[i] = 0

                err_tot=err_fs+err_ss_min+err_ss_max
                #print(err_ss_max)
                #print(err_ss_min)
                #print(err_fs)

        return  err_tot.sum()


    fitting_rule='monod'

    bound=False
    fitta=True
    load_fit=False#True
    plotta=False
    avanza=True
    modalita='Scrittura'#'Lettura'#

    t0_val=0
    n_test=200

    if modalita=='Scrittura':
        spt=[]
        spk_time_orig=[]
        spk_time=[]
        spk_time_orig_c=[]
        spk_time_c=[]
        file1 = open('dati_exp_'+neuron_nb+'.txt', 'r')
        Lines = file1.readlines()
        EL =np.float16(Lines[0])
        vrm =np.float16(Lines[1])/ -EL
        vth =np.float16(Lines[2])/ -EL
        Istm = np.int32(Lines[3].split(','))
        for i in range(len(Istm)):
            try:
                spk_time_orig.append(np.float64(Lines[4+i].split(',')))
            except:
                spk_time_orig.append([])

    #spk_time_orig.append([29,57,93,145,215,289,364,439, 514, 588 ,663, 738, 813 ])
    #spk_time_orig.append([30, 60, 100, 165, 246, 330, 413, 497, 580, 664, 748 ])
    #spk_time_orig.append([30, 63, 112, 194, 288, 382, 477, 572 ,667 ,762 ])
    #spk_time_orig.append([31, 67, 130, 235,344, 454, 564, 674 ,785])
    #spk_time_orig.append([32, 73,166, 298, 432, 566, 701 ])
    #spk_time_orig.append([34,83,267 ])
    #spk_time_orig.append([35,106])
    #spk_time_orig.append([37])
    #spk_time_orig.append([40])
        input_start_time=400
        dur_sign=400

        for i in range(len(spk_time_orig)):
            spk_time.append(np.array(spk_time_orig[i]) - input_start_time)

        spk_interval = []
        spk_interval_c = []
        spk_tr = []

        vol_tr = []


        for i in range(len(spk_time)):
            if spk_time[i].size > 0:
                spk_interval.append([spk_time[i][0]])
                corr_sp = i+1
                if spk_time[i].size > 1:
                    corr_sp2 = i+1
            else:
                spk_interval.append([])

        for j in range(len(spk_time_orig)):
            for i in range(len(spk_time[j]) - 1):
                spk_interval[j].append(spk_time[j][i + 1] - spk_time[j][i]-2)

        spk_time_orig_c=spk_time_orig
        corr_sp_c=corr_sp
        corr_sp2_c=corr_sp2
        Istm_c=Istm
        spk_time_c=spk_time



        varbound=np.array([[0,3000]]*6)

        if len(Istm)-corr_sp==0:
            varbound[2][0] = 300+0.0001  #valore minimo ith se il neurone spara per tutte le correnti
        else:
            varbound[2][0]=Istm[len(Istm)-corr_sp-1]+0.0001 
        varbound[2][1] =Istm[len(Istm)-corr_sp]-0.0001 #valore massimo ith (per la prima corrente che abbiamo spikes ith<300

        varbound[0][1] = 100000  # par_sc
        varbound[1][1] = 100000  # tao_m
        varbound[3][1] = 100  # Idep0 max
        if cm_fixed:
          varbound[4][0] = 189.79  # cm min
          varbound[4][1] = 189.79  # cm max
        varbound[5][1] = 1  # Idep0_ini max
        


        algorithm_param = {'max_num_iteration': 250,\
                           'population_size':200,\
                           'mutation_probability':0.3,\
                           'elit_ratio': 0.02,\
                           'crossover_probability': 0.3,\
                           'parents_portion': 0.2,\
                           'crossover_type':'uniform',\
                           'max_iteration_without_improv':None}

        model=ga(function=loss_func,dimension=6,variable_type='real',variable_boundaries=varbound,algorithm_parameters=algorithm_param,convergence_curve=False)
        print(model.param, neuron_nb)
        model.run()

        Vconvfact=-EL
        par_sc = model.best_variable[0]
        tao_m = model.best_variable[1]
        Ith = model.best_variable[2]
        Idep_ini_vr = model.best_variable[3]
        Cm=model.best_variable[4]
        Idep_ini=model.best_variable[5]

        t = sym.Symbol('t')
        delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')
        [V, Iadap, Idep] = load_v3()
        min_sc = (-Cm * EL) / tao_m + (2 / (1 + vth)) * (Ith + np.sqrt(Ith * (Ith - Cm * EL * (1 + vth) / tao_m)))
        sc = min_sc + par_sc
        k_adap_min = sc * (-EL / tao_m + Ith / (Cm * (1 + vth))) / (EL ** 2)
        k_adap_max = ((Cm * EL - sc * tao_m) ** 2) / (4 * Cm * (EL ** 2) * (tao_m ** 2))
        k_adap = k_adap_min + (k_adap_max - k_adap_min) * 0.01  #
        delta1 = -Cm * EL / (sc * tao_m)
        bet = Cm * (EL ** 2) * k_adap / (sc ** 2)
        psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)
        time_scale = 1 / (-sc / (Cm * EL))
        aux = Istm_c / sc
        alpha3 = aux.tolist()







        spt_00=np.zeros(len(Istm))
        for i in range(len(Istm)):
            aux=V.subs(alpha, alpha3[i]).subs(beta, bet).subs(delta, delta1).subs(t0,0).subs(V0,-1).subs(IaA0, 0).subs(IdA0, Idep_ini*(Istm[i]-Ith)*(Istm[i]>Ith)).subs(Psi,psi1)
            lam_x = sym.lambdify(t, aux, modules=['numpy'])
            x_vals = np.linspace(0, 1000, 10001)/ time_scale
            y_vals = lam_x(x_vals)
            aus = np.nonzero((y_vals > vth) * y_vals)
            if aus[0].size > 0:
                spt_00[i] = x_vals[aus[0][0]] * time_scale
            else:
                spt_00[i] = -1

        IaA_min = 0
        m = 1000


        time_sp=np.zeros([len(Istm),m])
        time_var =np.zeros([len(Istm),m-1])
        t_sp_tutti=[]
        t_sp_abs_tutti=[]
        Iada_tutti=[]
        t_var_tutti =[]
        time_soglia=[]
        spk_interval_modified = spk_interval

        for i in range(len(Istm)):
                IaA_max =  Idep_ini_vr + alpha3[i] / bet+(delta1/bet)*(1+vrm)

                IaA = np.linspace(IaA_min, IaA_max, m)
                t_sp = []
                t_sp_abs = []
                Iada=[]
                t_var = []

                for j in range(len(IaA)):
                    aux = V.subs(alpha, alpha3[i]).subs(beta, bet).subs(delta, delta1).subs(t0, 0).subs(V0, vrm).subs(IaA0,IaA[j]).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)

                    lam_x = sym.lambdify(t, aux, modules=['numpy'])
                    x_vals = np.linspace(0, 1000, 10001)/ time_scale
                    y_vals = lam_x(x_vals)
                    aus = np.nonzero((y_vals > vth) * y_vals)
                    if aus[0].size > 0:
                        time_sp[i,j] = x_vals[aus[0][0]] * time_scale
                        if j>0:
                            time_var[i,j-1] = (time_sp[i,j]-time_sp[i,j-1])/(IaA[j]-IaA[j-1])
                    else:
                        time_sp[i,j] = -1

                if len(spk_interval[len(Istm)-1-i])>1:
                    aux=spk_interval[len(Istm) - 1 - i][1]-(spt_00[i]-spk_interval[len(Istm) - 1 - i][0])
                    spk_interval_modified[len(Istm)-1-i][0]=spt_00[i]
                    spk_interval_modified[len(Istm) - 1 - i][1]=aux
                    t_sp.append(spt_00[i])
                    for j in range(len(spk_interval[len(Istm)-1-i])-1):
                        t_sp.append(time_sp[i,abs(time_sp[i, :] - spk_interval_modified[len(Istm)-1 - i][j+1]).argmin()])
                        Iada.append(IaA[abs(time_sp[i, :] - spk_interval_modified[len(Istm)-1 - i][j+1]).argmin()])
                        t_var.append(time_var[i,abs(IaA[0:len(time_var[len(Istm)-1-i])]-Iada[j]).argmin() - 1])

                        if j==0:
                            t_sp_abs.append(t_sp[j])
                        else:
                            t_sp_abs.append(t_sp_abs[j-1]+t_sp[j])

                    t_sp_tutti.append(t_sp)
                    t_sp_abs_tutti.append(t_sp_abs)
                    time_soglia.append(max(t_sp_abs)+(t_sp_abs[len(t_sp_abs)-1]-t_sp_abs[len(t_sp_abs)-2])/2)
                    Iada_tutti.append(Iada)
                    t_var_tutti.append(t_var)
                    if plotta:
                        plt.figure(170);
                        plt.scatter(t_sp_abs[0:len(t_sp_abs)], Iada)
                        plt.title('Iada_tot Idep'+str(Idep_ini_vr)+ ' tao ' + str(tao_m)+' Ith ' + str(Ith))
                        plt.figure(72);
                        plt.scatter(mat_to_iaa(np.linspace(0, len(time_sp[i]), len(time_sp[i]))), time_sp[i])
                        plt.title('time_su_iada_tot Idep' + str(Idep_ini_vr) + ' tao ' + str(tao_m) + ' Ith ' + str(Ith))
                        plt.xlabel('Iadap')
                        plt.ylabel('Time')
                    print("tempi spikes, Iadap")
                    print(t_sp_abs, Iada)
                else:
                    time_soglia.append(spt_00[i]/2)
        #break

    else:
        with open('neuron_' + neuron_nb + '_.pkl', 'rb') as f:  # Python 3: open(..., 'wb')
            [t_sp_abs_tutti, Iada_tutti, Istm_c, spk_time_c, t_sp_sim_c, Vconvfact, vth, vrm, bet, delta1, sc, time_scale,Idep_ini_vr, input_start_time, dur_sign, fitta, popt, corr_sp_c, corr_sp2_c, tao_m, Ith,time_soglia,Idep_ini,model] = pickle.load(f)
            corr_sp2=corr_sp2_c
            corr_sp=corr_sp_c
            Istm=Istm_c

            fitta = True
    xdata = np.linspace(0, int(input_start_time+dur_sign)+30, 100)

    popt_tutti=[]
    t_sp_sim_c=[]
        #a=[0,-250,-5,-5]
        #b=[1000,0,5,5]
        #prelast
        #a = [0, -1 / (Istm[len(Istm) - 1] / 1000), 0]
        #b = [200, 0, 2 / (Istm[len(Istm) - 1] / 1000)]

        #last
    if fitting_rule == 'exp_cum':
        a = [0, -1/(Istm[len(Istm)-1]/1000), 0,0]
        b=[500,0,2/(Istm[len(Istm)-1]/1000),1000]

        # monod
    if fitting_rule=='monod':
        a =  [0, 0, -200, 0]
        b =  [20000,10,200,100000]
        fun_loss_sel=monod_tot
        fun_sel = monod
    if bound:
        bnds = Bounds(np.array(a),np.array(b))

        #with open('best_conf.pkl','rb') as f:  # Python 3: open(..., 'rb')
        #    res,popt,Vconvfact, vtm, Vconvfact, vrm, popt, bet, delta1, alpha3, corr_sp2,sc, tao_m, Idep_ini_vr, input_start_time, dur_sign,is_Neuron_tr,time_soglia = pickle.load(f)


    if fitta:
        t_sp_sim_c=[]
        ada_sim_c=[]
        met='Nelder-Mead'#'slsqp'
        loss=[]
        if load_fit:
            with open('best_res_' + neuron_nb + '_cm.pkl', 'rb') as f:
                [best_res] = pickle.load(f)
        else:
            best_loss=np.inf
            for i in range(n_test):

                #aux=np.array(a)+np.random.rand(3)*(np.array(b)-np.array(a))
                aux = np.array(a) + np.random.rand(len(a)) * (np.array(b) - np.array(a))
                if bound:
                    res=minimize(fun_loss_sel,aux , method=met,bounds=bnds,options={'maxiter': 50000, 'disp': True})
                else:
                    res = minimize(fun_loss_sel, aux, method=met, options={'maxiter': 50000, 'disp': True})
                loss.append(res.fun)
                if(res.fun<best_loss):
                    best_loss=res.fun
                    best_res=res
            if cm_fixed:
                with open('best_res_' + neuron_nb + '_.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
                    pickle.dump([best_res], f)
            else:
                with open('best_res_' + neuron_nb + '_cm.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
                        pickle.dump([best_res], f)
        #popt = np.array([0., 0.])



        popt = np.array([0., 0.,0.,0.])
        for i in range(corr_sp_c):#
            min_fit = np.array([0.0, 0.0])
            max_fit = np.array([1000.0, 2.0])
            guess = np.array([1, 0.1])
                    #print(Istm_c[len(Istm_c)-1-i])

                    #popt[0] = best_res.x[0] + best_res.x[1] * ((6-corr_sp+i)*0.2)
                    #popt[1] = best_res.x[2] + best_res.x[3] * ((6-corr_sp+ i)*0.2)
            if fitting_rule=='monod':
                popt[0] = best_res.x[0]
                popt[1] = best_res.x[1] * Istm_c[len(Istm_c)-corr_sp_c+i]/1000
                popt[2] = best_res.x[2]
                popt[3] = best_res.x[3]


            aux = Istm_c / sc
            alpha=aux.tolist()
            is_Neuron_tr = True
            [t_aux,ada_aux,time,voltage]=plot_tr_from_fit_neuron_dep_ini2(Vconvfact, vth, vrm, fun_sel , popt, bet, delta1, alpha[len(alpha)-corr_sp_c+i],sc, time_scale, Idep_ini_vr, input_start_time, dur_sign,is_Neuron_tr,time_soglia[len(alpha)-corr_sp_c+i],Idep_ini,Ith)
            t_sp_sim_c.append(t_aux)
            ada_sim_c.append(ada_aux)

            popt_tutti.append(popt)
            if len(spk_time_c[len(Istm)-1-i])>1:
                if plotta:
                    plt.figure(170);

                #plt.plot(xdata, exp_cum(xdata, *popt), label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
                if popt.size == 2:
                    if plotta:
                        plt.plot(xdata, fun_sel(xdata, *popt), label='fit: a=%5.3f, b=%5.3f' % tuple(popt))

                if popt.size == 3:
                    if plotta:
                        plt.plot(xdata, fun_sel(xdata, *popt), label='fit: a=%5.3f, b=%5.3f,c=%5.3f' % tuple(popt))

                if popt.size==4:
                    if plotta:
                        plt.plot(xdata, fun_sel(xdata, *popt), label='fit: a=%5.3f, b=%5.3f,c=%5.3f,alpha=%5.3f' % tuple(popt))
                if plotta:
                    plt.legend(bbox_to_anchor=(0.2, 0.2), loc='upper left', borderaxespad=0.)
                    plt.xlabel('Time')
                    plt.ylabel('Iadap')

        #if cm_fixed:
        #    MyFile = open('Iadap_tsp' + neuron_nb + '.txt', 'w')
        #else:
        #    MyFile = open('Iadap_tsp'+neuron_nb+'_cm.txt', 'w')

#        for element in ada_sim_c:
#            MyFile.write(str(element))
#            MyFile.write('\n')
#        MyFile.close()
    else:
        popt = np.array([0., 0.,0.,0.])
        for i in range(corr_sp2_c):#
            min_fit = np.array([0.0, 0.0])
            max_fit = np.array([1000.0, 2.0])
            guess = np.array([1, 0.1])

            aux = Istm_c / sc
            alpha=aux.tolist()
            [t_aux,time,voltage]=plot_tr_v3_vect3_dep_ini2(Vconvfact, vth, vrm, Iada_tutti[i], bet, delta1, alpha[len(alpha)-corr_sp2_c+i],sc, time_scale, Idep_ini_vr, input_start_time, dur_sign,Idep_ini,Ith)
            t_sp_sim_c.append(t_aux)
            popt_tutti.append(popt)

    for i in range(corr_sp_c):
        if plotta:
            plt.figure(200);
            plt.scatter(spk_time_c[i], Istm_c[len(Istm_c)-1-i] * np.ones(len(spk_time_c[i])), marker='|',color='r');
            plt.scatter(t_sp_sim_c[len(t_sp_sim_c)-1-i],Istm_c[len(Istm_c)-1-i] * np.ones(len(t_sp_sim_c[len(t_sp_sim_c)-1-i])), marker='|', color='b');
            plt.ylabel('Current')
            plt.xlabel('spike times')

            plt.figure(201)
            plt.scatter(Istm_c[len(Istm_c)-1-i],len(spk_time_c[i]),marker='*',color='r')
            plt.scatter(Istm_c[len(Istm_c)-1-i],len(t_sp_sim_c[len(t_sp_sim_c)-1-i]),marker='*',color='b')
            plt.ylabel('number of spikes')
            plt.xlabel('Current')
    #with open('neuron_' + neuron_nb + '_.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
    #    pickle.dump([t_sp_abs_tutti, Iada_tutti, Istm_c, spk_time_c, t_sp_sim_c, Vconvfact, vth, vrm, bet, delta1,sc, time_scale, Idep_ini_vr, input_start_time, dur_sign, fitta, popt, corr_sp_c, corr_sp2_c, tao_m, Ith], f)
    if modalita=='Scrittura':
        if cm_fixed:
            with open('neuron_' + neuron_nb + '_idep.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
                pickle.dump([t_sp_abs_tutti, Iada_tutti, Istm_c, spk_time_c, t_sp_sim_c, Vconvfact, vth, vrm, bet, delta1,sc, time_scale, Idep_ini_vr, input_start_time, dur_sign, fitta, popt, corr_sp_c, corr_sp2_c, tao_m, Ith,time_soglia,Idep_ini], f)
        else:
            with open('neuron_' + neuron_nb + '_cm_idep.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
                pickle.dump([t_sp_abs_tutti, Iada_tutti, Istm_c, spk_time_c, t_sp_sim_c, Vconvfact, vth, vrm, bet, delta1,sc, time_scale, Idep_ini_vr, input_start_time, dur_sign, fitta, popt, corr_sp_c, corr_sp2_c, tao_m, Ith,time_soglia,Idep_ini], f)

        if cm_fixed:
            file_info = open(neuron_nb + "_info.txt", mode="w", encoding="utf-8")
        else:
            file_info = open(neuron_nb + "_info_cm.txt", mode="w", encoding="utf-8")
        file_info.write('Cm='+str(Cm) + '\n' +'Ith='+str(Ith) + '\n' +'tao='+ str(tao_m) + '\n' +'sc='+ str(sc)+ '\n' +'alpha='+ str(alpha) + '\n'+ 'bet='+str(bet)+ '\n'+'delta1='+str(delta1)+ '\n'+'Idep_ini='+str(Idep_ini)+'\n'+'Idep_ini_vr='+str(Idep_ini_vr)+ '\n'+'psi='+str(psi1)+ '\n'+'time scale='+str(time_scale)+'\n'+'A='+str(best_res.x[0])+'\n'+'B='+str(best_res.x[1])+'\n'+'C='+str(best_res.x[2])+'\n'+'alpha='+str(best_res.x[3]))
        file_info.close()
    print("Ith,tao_m,sc,alpha3,bet,delta1,Idep_ini_vr,psi1,time scale")
    print(Ith,tao_m,sc,alpha,bet,delta1,Idep_ini_vr,psi1,time_scale)

    return Ith
