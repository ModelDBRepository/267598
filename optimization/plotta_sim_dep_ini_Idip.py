def plot_tr_v3_vect3_dep_ini(Vconvfact,vtm,vrm,a_inc,bet,delta1,alpha3,sc, tao_m, Idep_ini_vr, st_sign, dur_sign,Idep_ini,Ith):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    plotta=False
    n_sp=np.size(a_inc)

    corr = round(alpha3 * sc) / 1000

    tim_aux=np.linspace(0, 1000, 10000)

    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
    #cell_num = '95810005'

    t = sym.Symbol('t')
    delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load_v3()
    psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)

    t0_val=0
    #vtm=np.mean(vol_tr)/66.35
    #vrm=-1.0
    #d_in=st_point_dep1[r]
    #a_in=st_point_ada1[r]
    #aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0,d_in)
    #aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t0_val).subs(V0, -1).subs(IaA0,a_in).subs(IdA0,d_in).subs(Psi, psi1)
    if plotta:
        plt.figure()
    time=[]
    voltage=[]
    t_next=0
    adap_final=[]
    for i in range(n_sp+1):
        if i>0:
            t_init=t_next+2/tao_m
        else:
            t_init=t_next


        #aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if i==0:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0, Idep_ini*(corr*1000-Ith)*(corr*1000>Ith)).subs(Psi, psi1)
        else:
            print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            print('***** ada inc ****' + str(i) + ' *****spike****')
            print(a_inc[i-1])
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, a_inc[i-1]).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)

        #print('ada comp')
        #print( ada+a_inc)
        #print('dep comp')
        #print(dep + d_inc)

#        print('ada')
#        print(ada_nc)
#        print('dep')
#        print(dep_nc)
        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / tao_m)>t_init]/ tao_m
        y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:



            t_next=x_vals[y_vals>vtm][0]
            t_out.append(t_next*tao_m)

            if i > 0:
                adap_final.append(Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0,a_inc[i - 1]).subs(IdA0, Idep_ini_vr).subs(Psi, psi1).subs(t, t_next))
            else:
                adap_final.append(Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0, 0).subs(IdA0, Idep_ini*(corr*1000-Ith)*(corr*1000>Ith)).subs(Psi, psi1).subs(t, t_next))


#            ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc[i]).subs(IdA0, dep+d_inc[i]).subs(Psi, psi1).subs(t, t_next)
#            dep = Idep.subs(beta, bet).subs(IdA0, dep+d_inc[i]).subs(t, t_next).subs(t0, t_init)
            #ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
            #dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
#        auxa = (IaA - a_in) ** 2
#        auxd = (IdA - d_in) ** 2
#        ada_nc = Iadap_mat[4, auxa.argmin(), auxd.argmin()]
#        dep_nc = Idep_mat[4, auxa.argmin(), auxd.argmin()]
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            if plotta:
                plt.plot(st_sign + x_vals * tao_m, y_vals * Vconvfact)
            time.append(st_sign + x_vals * tao_m)
            voltage.append(y_vals * Vconvfact)
        else:
            print('aoa')
            print(n_sp)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)
            t_next = dur_sign/tao_m
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals
            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            if plotta:
                plt.plot(st_sign + x_vals * tao_m, y_vals * Vconvfact)
            time.append(st_sign + x_vals * tao_m)
            voltage.append(y_vals * Vconvfact)

            i=n_sp+1
        t_init = t_next + 2 / tao_m
    #plt.plot(tim, vol)
    #plt.title('trace ('+ str(sc*alpha3)+'nA)')
    return t_out,time,voltage,adap_final


def plot_tr_from_fit_neuron_rette_dep_ini(Vconvfact, vtm, vrm, funzione,popt, bet, delta1,alpha3, sc,tao,Idep_ini_vr,st_sign,dur_sign,lin_func_inf,vinc_inf,lin_func_sup,vinc_sup,Idep_ini,Ith):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3
    import sympy as sym
    I = sym.Symbol('I')
    plotta=False
    #n_sp=np.size(a_inc)

    corr=round(alpha3*sc)/1000
    if corr*1000<=vinc_inf:
        dur_sign=min(dur_sign,lin_func_inf.subs(I,corr*1000))

    if corr*1000>=vinc_sup:
        dur_sign=min(dur_sign,lin_func_sup.subs(I,corr*1000))
    tim = []


    tim_aux = np.linspace(0, 1000, 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    #cell_num = '95810005'

    t = sym.Symbol('t')
    delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load_v3()
    psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)

    t_init=0
    if plotta:
        plt.figure()
    i=0
    ada_vec=[]
    time=[]
    voltage=[]
    #while i<20:
    while t_init * tao < dur_sign :
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada  ****' + str(i) + ' *****spike****')
        if i==0:
            print(funzione(t_init * tao, *popt))
        else:
            new_ada=funzione(t_next*tao,*popt)
            ada_vec.append(new_ada)
            print(new_ada)
            print(funzione(t_next*tao,*popt))

        print('t_init')
        print(t_init)
        print(t_init * tao)


        # aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if t_init * tao== dur_sign:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0,adaf).subs(IdA0, depf).subs(Psi, psi1)

        else:
            if i==0:
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0,Idep_ini*(corr*1000-Ith)*(corr*1000>Ith)).subs(Psi, psi1)

            else:
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, new_ada).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / tao) > t_init] / tao
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
            if i == 0:
                adaf = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, 0).subs(IdA0, Idep_ini*(corr*1000-Ith)*(corr*1000>Ith)).subs(Psi, psi1).subs(t, t_next)
                depf = Idep.subs(beta, bet).subs(IdA0, Idep_ini*(corr*1000-Ith)*(corr*1000>Ith)).subs(t0, t_init).subs(t, t_next)
            else:
                adaf = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0,vrm).subs(IaA0, new_ada).subs(IdA0, Idep_ini_vr).subs(Psi, psi1).subs(t, t_next)
                depf = Idep.subs(beta, bet).subs(IdA0, Idep_ini_vr).subs(t0, t_init).subs(t, t_next)
            print('.................adap....')
            print(adaf)
            print('..................dep....')
            print(depf)
            if (t_next*tao<dur_sign):
                t_out.append(t_next*tao)
                print("t_next .....")
                print(t_next)
                print(t_next*tao)
            else:
                t_next=dur_sign/tao

            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals

            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            time.append(st_sign + x_vals * tao)
            voltage.append(y_vals * Vconvfact)
            if plotta:
                plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact)
            #st_sign=st_sign#+t_next*15.58
            print("st_sign .....")
            print(st_sign)
        else:
            print('aoa')
            print(i)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)
            t_next = dur_sign/tao
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals
            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            try:
                time.append(st_sign + x_vals * tao)
                voltage.append(y_vals * Vconvfact)
                if plotta:
                    plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact)
            except:
                time.append(st_sign + x_vals * tao)
                voltage.append(y_vals * Vconvfact*np.ones(x_vals.size))
                if plotta:
                    plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact*np.ones(x_vals.size))

        i = i + 1
        print('t_next')
        print(t_next)
        if t_next==dur_sign:
            adaf=Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0,new_ada ).subs(IdA0, Idep_ini_vr).subs(Psi, psi1).subs(t,dur_sign)
            depf=Idep.subs(beta, bet).subs(IdA0, Idep_ini_vr).subs(t0, t_init).subs(t,dur_sign)
            t_init=t_next
            vrm=y_vals[len(y_vals)-1]
            alpha3=0

        else:
            t_init = t_next + 2 / tao
    try:
        if plotta:
            plt.plot(tim, vol)
    except:
        print('no trace')
    if plotta:
        plt.title('trace (' + str(alpha3*sc/1000) + 'nA)')
    return t_out,ada_vec,time,voltage
