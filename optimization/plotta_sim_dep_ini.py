def plot_tr(vol_tr,st_point_dep1,st_point_ada1,dep_vec1_2,ada_vec1_2,bet,gam,alpha3,n_sp,r,c):
    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load

    import pyabf
    cell_num = '95810005'
    abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
    abf.setSweep(15)
    vol_base = -62.86698717948718
    abf.sweepY = abf.sweepY + vol_base

    t = sym.Symbol('t')
    alpha,beta,gamma,IaA0,IdA0,t0,V0= sym.symbols('alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load()
    #V=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*((-2)*sym.exp(1)**(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta*((-1)+gamma)+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(1/2)*((-1)+alpha+beta)*((-1)+beta+gamma**2)+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)*V0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+alpha*((-1)+(1+(-1)*beta)**(1/2)+beta)*((-1)+beta+gamma**2)+(-1)*((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+((-1)+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(-1)*((-1)+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2)))+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(alpha*(1+(1+(-1)*beta)**(1/2)+(-1)*beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+(-1)*(1+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(1+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2))))
    #Iadap=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*(2*sym.exp(1) **(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*alpha*(1+(-1)*beta)**(1/2)*((-1)+beta+gamma**2)+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*alpha*(1+(-1)*beta)**(1/2)+(-2)*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+IdA0*beta*gamma+(-1)*IdA0*beta**2*gamma+(-1)*gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+beta*gamma**2+(-1)*IaA0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2)))+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)+(-1)*alpha*(1+(-1)*beta)**(1/2)+2*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+(-1)*beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+(-1)*IdA0*beta*gamma+IdA0*beta**2*gamma+gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+(-1)*beta*gamma**2+(-1)*IaA0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2))));
    #Idep=sym.exp(1)**(((-1)*t+t0)*gamma)*IdA0;





    t0_val=0
    vtm=np.mean(vol_tr)/66.35
    vrm=-1.0
    d_in=st_point_dep1[r]
    a_in=st_point_ada1[r]
    aux = V.subs(alpha, alpha3[1]).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, vrm).subs(IaA0, a_in).subs(IdA0,d_in)
    print('***** dep init ****')
    print(d_in)
    print('***** ada init ****')
    print(a_in)

    lam_x = sym.lambdify(t, aux, modules=['numpy'])
    x_vals = abf.sweepX * 1000 / 15.58
    y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
    t_next=x_vals[y_vals>vtm][0]

    ada = Iadap.subs(alpha, alpha3[1]).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, vrm).subs(IaA0, a_in).subs(IdA0, d_in).subs(t,t_next)
    dep = Idep.subs(gamma, gam).subs(t0, t0_val).subs(IdA0, d_in).subs(t,t_next)


    #with open('spt_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    spt = pickle.load(f)

    #with open('Iadap_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    Iadap_mat = pickle.load(f)

    #with open('Idep_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + 'gam_' + str(gam) + '_.pkl', "rb") as f:
    #    Idep_mat = pickle.load(f)

    #auxa = (IaA - a_in) ** 2
    #auxd = (IdA - d_in) ** 2
    #ada_nc=Iadap_mat[4,auxa.argmin(),auxd.argmin()]
    #dep_nc=Idep_mat[4, auxa.argmin(),auxd.argmin()]
    x_sector = np.logical_and(t0_val < x_vals, x_vals < t_next) * x_vals


    x_vals = x_sector[np.nonzero(x_sector)]
    y_vals = lam_x(x_vals)
    plt.figure()
    plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


    vrm=-0.7349960338647149
    d_inc=dep_vec1_2[r,c]
    a_inc=ada_vec1_2[r,c]
    for i in range(n_sp):
        t_init=t_next+2/15.58
        print('******dep ****'+str(i)+' *****spike****')
        print(dep + d_inc)
        print('***** ada ****'+str(i)+' *****spike****')
        print(ada + a_inc)


        aux = V.subs(alpha, alpha3[1]).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        #print('ada comp')
        #print( ada+a_inc)
        #print('dep comp')
        #print(dep + d_inc)

#        print('ada')
#        print(ada_nc)
#        print('dep')
#        print(dep_nc)
        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = abf.sweepX[(abf.sweepX * 1000 / 15.58)>t_init]* 1000 / 15.58
        y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        t_init = t_next + 2 / 15.58

        if len(np.nonzero(y_vals > vtm)[0]) > 0:
            t_next = x_vals[y_vals > vtm][0]


            t_next=x_vals[y_vals>vtm][0]
            ada = Iadap.subs(alpha, alpha3[1]).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
            dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
#        auxa = (IaA - a_in) ** 2
#        auxd = (IdA - d_in) ** 2
#        ada_nc = Iadap_mat[4, auxa.argmin(), auxd.argmin()]
#        dep_nc = Idep_mat[4, auxa.argmin(), auxd.argmin()]
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            #plt.figure()
            plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        else:
            i=n_sp+1

def plot_tr2(vtm,a_in,d_in,a_inc,d_inc,bet,gam,n_sp,cell_num,alpha3):
    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/0.252587)+10
    #cell_num = '95810005'
    abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
    abf.setSweep(tr)
    vol_base = -62.86698717948718
    abf.sweepY = abf.sweepY + vol_base

    t = sym.Symbol('t')
    alpha,beta,gamma,IaA0,IdA0,t0,V0= sym.symbols('alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load()
    #V=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*((-2)*sym.exp(1)**(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta*((-1)+gamma)+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(1/2)*((-1)+alpha+beta)*((-1)+beta+gamma**2)+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)*V0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+alpha*((-1)+(1+(-1)*beta)**(1/2)+beta)*((-1)+beta+gamma**2)+(-1)*((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+((-1)+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(-1)*((-1)+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2)))+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(alpha*(1+(1+(-1)*beta)**(1/2)+(-1)*beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+(-1)*(1+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(1+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2))))
    #Iadap=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*(2*sym.exp(1) **(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*alpha*(1+(-1)*beta)**(1/2)*((-1)+beta+gamma**2)+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*alpha*(1+(-1)*beta)**(1/2)+(-2)*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+IdA0*beta*gamma+(-1)*IdA0*beta**2*gamma+(-1)*gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+beta*gamma**2+(-1)*IaA0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2)))+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)+(-1)*alpha*(1+(-1)*beta)**(1/2)+2*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+(-1)*beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+(-1)*IdA0*beta*gamma+IdA0*beta**2*gamma+gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+(-1)*beta*gamma**2+(-1)*IaA0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2))));
    #Idep=sym.exp(1)**(((-1)*t+t0)*gamma)*IdA0;





    t0_val=0
    #vtm=np.mean(vol_tr)/66.35
    vrm=-1.0
    #d_in=st_point_dep1[r]
    #a_in=st_point_ada1[r]
    aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, vrm).subs(IaA0, a_in).subs(IdA0,d_in)
    print('***** dep init ****')
    print(d_in)
    print('***** ada init ****')
    print(a_in)

    lam_x = sym.lambdify(t, aux, modules=['numpy'])
    x_vals = abf.sweepX * 1000 / 15.58
    y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
    t_next=x_vals[y_vals>vtm][0]
    t_out.append(t_next)

    ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, vrm).subs(IaA0, a_in).subs(IdA0, d_in).subs(t,t_next)
    dep = Idep.subs(gamma, gam).subs(t0, t0_val).subs(IdA0, d_in).subs(t,t_next)


    #with open('spt_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    spt = pickle.load(f)

    #with open('Iadap_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    Iadap_mat = pickle.load(f)

    #with open('Idep_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + 'gam_' + str(gam) + '_.pkl', "rb") as f:
    #    Idep_mat = pickle.load(f)

    #auxa = (IaA - a_in) ** 2
    #auxd = (IdA - d_in) ** 2
    #ada_nc=Iadap_mat[4,auxa.argmin(),auxd.argmin()]
    #dep_nc=Idep_mat[4, auxa.argmin(),auxd.argmin()]
    x_sector = np.logical_and(t0_val < x_vals, x_vals < t_next) * x_vals


    x_vals = x_sector[np.nonzero(x_sector)]
    y_vals = lam_x(x_vals)
    plt.figure()
    plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


    vrm=-0.7349960338647149
    #d_inc=dep_vec1_2[r,c]
    #a_inc=ada_vec1_2[r,c]
    for i in range(n_sp):
        t_init=t_next+2/15.58
        print('******dep ****'+str(i)+' *****spike****')
        print(dep + d_inc)
        print('***** ada ****'+str(i)+' *****spike****')
        print(ada + a_inc)


        aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        #print('ada comp')
        #print( ada+a_inc)
        #print('dep comp')
        #print(dep + d_inc)

#        print('ada')
#        print(ada_nc)
#        print('dep')
#        print(dep_nc)
        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = abf.sweepX[(abf.sweepX * 1000 / 15.58)>t_init]* 1000 / 15.58
        y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        t_init = t_next + 2 / 15.58

        if len(np.nonzero(y_vals > vtm)[0]) > 0:



            t_next=x_vals[y_vals>vtm][0]
            t_out.append(t_next)
            ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
            dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
#        auxa = (IaA - a_in) ** 2
#        auxd = (IdA - d_in) ** 2
#        ada_nc = Iadap_mat[4, auxa.argmin(), auxd.argmin()]
#        dep_nc = Idep_mat[4, auxa.argmin(), auxd.argmin()]
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            #plt.figure()
            plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        else:
            i=n_sp+1

    plt.plot(abf.sweepX * 1000, abf.sweepY)
    return t_out


def plot_tr3(vtm,a_in,d_in,a_inc,d_inc,bet,gam,n_sp,cell_num,alpha3):
    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    tr=min(int(alpha3/0.252587)+10,15)
    #cell_num = '95810005'
    abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
    abf.setSweep(tr)
    vol_base = -62.86698717948718
    abf.sweepY = abf.sweepY + vol_base

    t = sym.Symbol('t')
    alpha,beta,gamma,IaA0,IdA0,t0,V0= sym.symbols('alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load()
    #V=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*((-2)*sym.exp(1)**(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta*((-1)+gamma)+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(1/2)*((-1)+alpha+beta)*((-1)+beta+gamma**2)+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)*V0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+alpha*((-1)+(1+(-1)*beta)**(1/2)+beta)*((-1)+beta+gamma**2)+(-1)*((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+((-1)+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(-1)*((-1)+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2)))+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(alpha*(1+(1+(-1)*beta)**(1/2)+(-1)*beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+(-1)*(1+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(1+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2))))
    #Iadap=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*(2*sym.exp(1) **(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*alpha*(1+(-1)*beta)**(1/2)*((-1)+beta+gamma**2)+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*alpha*(1+(-1)*beta)**(1/2)+(-2)*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+IdA0*beta*gamma+(-1)*IdA0*beta**2*gamma+(-1)*gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+beta*gamma**2+(-1)*IaA0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2)))+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)+(-1)*alpha*(1+(-1)*beta)**(1/2)+2*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+(-1)*beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+(-1)*IdA0*beta*gamma+IdA0*beta**2*gamma+gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+(-1)*beta*gamma**2+(-1)*IaA0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2))));
    #Idep=sym.exp(1)**(((-1)*t+t0)*gamma)*IdA0;





    t0_val=0
    #vtm=np.mean(vol_tr)/66.35
    vrm=-1.0
    #d_in=st_point_dep1[r]
    #a_in=st_point_ada1[r]
    aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, vrm).subs(IaA0, a_in).subs(IdA0,d_in)
    print('***** dep init ****')
    print(d_in)
    print('***** ada init ****')
    print(a_in)

    lam_x = sym.lambdify(t, aux, modules=['numpy'])
    x_vals = abf.sweepX * 1000 / 15.58
    y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
    t_out=[]
    if len(np.nonzero(y_vals > vtm)[0]) > 0:
        t_next=x_vals[y_vals>vtm][0]
        t_out.append(t_next)
        ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, vrm).subs(IaA0, a_in).subs(IdA0, d_in).subs(t,t_next)
        dep = Idep.subs(gamma, gam).subs(t0, t0_val).subs(IdA0, d_in).subs(t,t_next)


        #with open('spt_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
        #    spt = pickle.load(f)

        #with open('Iadap_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
        #    Iadap_mat = pickle.load(f)

        #with open('Idep_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + 'gam_' + str(gam) + '_.pkl', "rb") as f:
        #    Idep_mat = pickle.load(f)

        #auxa = (IaA - a_in) ** 2
        #auxd = (IdA - d_in) ** 2
        #ada_nc=Iadap_mat[4,auxa.argmin(),auxd.argmin()]
        #dep_nc=Idep_mat[4, auxa.argmin(),auxd.argmin()]
        x_sector = np.logical_and(t0_val < x_vals, x_vals < t_next) * x_vals


        x_vals = x_sector[np.nonzero(x_sector)]
        y_vals = lam_x(x_vals)
        plt.figure()
        plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        vrm=-0.7349960338647149
        #d_inc=dep_vec1_2[r,c]
        #a_inc=ada_vec1_2[r,c]
        for i in range(n_sp):
            t_init=t_next+2/15.58
            print('******dep ****'+str(i)+' *****spike****')
            print(dep + d_inc)
            print('***** ada ****'+str(i)+' *****spike****')
            print(ada + a_inc)


            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
            #print('ada comp')
            #print( ada+a_inc)
            #print('dep comp')
            #print(dep + d_inc)

    #        print('ada')
    #        print(ada_nc)
    #        print('dep')
    #        print(dep_nc)
            lam_x = sym.lambdify(t, aux, modules=['numpy'])
            x_vals = abf.sweepX[(abf.sweepX * 1000 / 15.58)>t_init]* 1000 / 15.58
            y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
            t_init = t_next + 2 / 15.58

            if len(np.nonzero(y_vals > vtm)[0]) > 0:

                t_next=x_vals[y_vals>vtm][0]
                t_out.append(t_next)
                ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
                dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
    #        auxa = (IaA - a_in) ** 2
    #        auxd = (IdA - d_in) ** 2
    #        ada_nc = Iadap_mat[4, auxa.argmin(), auxd.argmin()]
    #        dep_nc = Idep_mat[4, auxa.argmin(), auxd.argmin()]
                x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


                x_vals = x_sector[np.nonzero(x_sector)]
                y_vals = lam_x(x_vals)
                #plt.figure()
                plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
            else:
                i=n_sp+1

        plt.plot(abf.sweepX * 1000, abf.sweepY)
    else:
        plt.figure()
        plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        plt.plot(abf.sweepX * 1000, abf.sweepY)
    return t_out




def plot_trnew(Vconvfact,vtm,vrm,a_in,d_in,a_inc,d_inc,bet,gam,n_sp,cell_num,alpha3):
    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/0.252587)+10
    #cell_num = '95810005'
    abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
    abf.setSweep(tr)
    vol_base = -62.86698717948718
    abf.sweepY = abf.sweepY + vol_base

    t = sym.Symbol('t')
    alpha,beta,gamma,IaA0,IdA0,t0,V0= sym.symbols('alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load()
    #V=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*((-2)*sym.exp(1)**(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta*((-1)+gamma)+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(1/2)*((-1)+alpha+beta)*((-1)+beta+gamma**2)+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)*V0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+alpha*((-1)+(1+(-1)*beta)**(1/2)+beta)*((-1)+beta+gamma**2)+(-1)*((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+((-1)+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(-1)*((-1)+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2)))+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(alpha*(1+(1+(-1)*beta)**(1/2)+(-1)*beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+(-1)*(1+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(1+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2))))
    #Iadap=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*(2*sym.exp(1) **(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*alpha*(1+(-1)*beta)**(1/2)*((-1)+beta+gamma**2)+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*alpha*(1+(-1)*beta)**(1/2)+(-2)*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+IdA0*beta*gamma+(-1)*IdA0*beta**2*gamma+(-1)*gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+beta*gamma**2+(-1)*IaA0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2)))+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)+(-1)*alpha*(1+(-1)*beta)**(1/2)+2*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+(-1)*beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+(-1)*IdA0*beta*gamma+IdA0*beta**2*gamma+gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+(-1)*beta*gamma**2+(-1)*IaA0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2))));
    #Idep=sym.exp(1)**(((-1)*t+t0)*gamma)*IdA0;





    t0_val=0
    #vtm=np.mean(vol_tr)/66.35
    #vrm=-1.0
    #d_in=st_point_dep1[r]
    #a_in=st_point_ada1[r]
    aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0,d_in)
    print('***** dep init ****')
    print(d_in)
    print('***** ada init ****')
    print(a_in)

    lam_x = sym.lambdify(t, aux, modules=['numpy'])
    x_vals = abf.sweepX * 1000 / 15.58
    y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
    t_next=x_vals[y_vals>vtm][0]
    t_out.append(t_next)

    ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0, d_in).subs(t,t_next)
    dep = Idep.subs(gamma, gam).subs(t0, t0_val).subs(IdA0, d_in).subs(t,t_next)


    #with open('spt_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    spt = pickle.load(f)

    #with open('Iadap_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    Iadap_mat = pickle.load(f)

    #with open('Idep_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + 'gam_' + str(gam) + '_.pkl', "rb") as f:
    #    Idep_mat = pickle.load(f)

    #auxa = (IaA - a_in) ** 2
    #auxd = (IdA - d_in) ** 2
    #ada_nc=Iadap_mat[4,auxa.argmin(),auxd.argmin()]
    #dep_nc=Idep_mat[4, auxa.argmin(),auxd.argmin()]
    x_sector = np.logical_and(t0_val < x_vals, x_vals < t_next) * x_vals


    x_vals = x_sector[np.nonzero(x_sector)]
    y_vals = lam_x(x_vals)
    plt.figure()
    plt.plot(31.1 + x_vals * 15.58, y_vals * Vconvfact)


    #vrm=-0.7349960338647149
    #d_inc=dep_vec1_2[r,c]
    #a_inc=ada_vec1_2[r,c]
    for i in range(n_sp):
        t_init=t_next+2/15.58
        print('******dep ****'+str(i)+' *****spike****')
        print(dep + d_inc)
        print('***** ada ****'+str(i)+' *****spike****')
        print(ada + a_inc)


        aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        #print('ada comp')
        #print( ada+a_inc)
        #print('dep comp')
        #print(dep + d_inc)

#        print('ada')
#        print(ada_nc)
#        print('dep')
#        print(dep_nc)
        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = abf.sweepX[(abf.sweepX * 1000 / 15.58)>t_init]* 1000 / 15.58
        y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        t_init = t_next + 2 / 15.58

        if len(np.nonzero(y_vals > vtm)[0]) > 0:



            t_next=x_vals[y_vals>vtm][0]
            t_out.append(t_next)
            ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
            dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
#        auxa = (IaA - a_in) ** 2
#        auxd = (IdA - d_in) ** 2
#        ada_nc = Iadap_mat[4, auxa.argmin(), auxd.argmin()]
#        dep_nc = Idep_mat[4, auxa.argmin(), auxd.argmin()]
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            #plt.figure()
            plt.plot(31.1 + x_vals * 15.58, y_vals * Vconvfact)
        else:
            print('aoa')
            print(n_sp)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)

            i=n_sp+1

    #plt.plot(abf.sweepX * 1000, abf.sweepY)
    return t_out



def plot_trnew2(Vconvfact,vtm,vrm,a_in,d_in,a_inc,d_inc,bet,gam,n_sp,cell_num,alpha3):
    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load

    st_sign = 531.1
    corr=int(alpha3/0.252587)*0.2
    f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

    aux=f.readline()
    tim=[]
    vol=[]
    trovato = False
    for aux in f:
        tim.append(np.float64(aux.split('\t')[0]))
        vol.append(np.float64(aux.split('\t')[1]))
    tim = np.array(tim)
    vol = np.array(vol)

    tim_aux=np.linspace(tim.min(), tim.max(), 10000)

    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/0.252587)+10
    #cell_num = '95810005'


    t = sym.Symbol('t')
    alpha,beta,gamma,IaA0,IdA0,t0,V0= sym.symbols('alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load()
    #V=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*((-2)*sym.exp(1)**(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta*((-1)+gamma)+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(1/2)*((-1)+alpha+beta)*((-1)+beta+gamma**2)+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)*V0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+alpha*((-1)+(1+(-1)*beta)**(1/2)+beta)*((-1)+beta+gamma**2)+(-1)*((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+((-1)+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(-1)*((-1)+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2)))+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(alpha*(1+(1+(-1)*beta)**(1/2)+(-1)*beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+((-1)+beta)*((IaA0+(-1)*IdA0)*beta**2+(-1)*(1+(1+(-1)*beta)**(1/2))*beta*((-1)+IdA0*((-1)+gamma))+(1+(1+(-1)*beta)**(1/2))*((-1)+gamma**2)+IaA0*beta*((-1)+gamma**2))))
    #Iadap=(1/2)*sym.exp(1)**((-1)*t0*(1+(-1)*beta)**(1/2)+(-1)*t*((1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*beta)**(-3/2)*((-1)+beta+gamma**2)**(-1)*(2*sym.exp(1) **(t*(1+(-1)*beta)**(1/2)+t0*((1+(-1)*beta)**(1/2)+gamma))*IdA0*(1+(-1)*beta)**(3/2)*beta+(-2)*sym.exp(1)**(t0*(1+(-1)*beta)**(1/2)+t*((1+(-1)*beta)**(1/2)+gamma))*alpha*(1+(-1)*beta)**(1/2)*((-1)+beta+gamma**2)+sym.exp(1)**(t*(2*(1+(-1)*beta)**(1/2)+gamma))*(1+(-1)*alpha*(1+(-1)*beta)**(1/2)+(-2)*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+IdA0*beta*gamma+(-1)*IdA0*beta**2*gamma+(-1)*gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+beta*gamma**2+(-1)*IaA0*((-1)+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+(-1)*V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2)))+sym.exp(1)**(2*t0*(1+(-1)*beta)**(1/2)+t*gamma)*((-1)+(-1)*alpha*(1+(-1)*beta)**(1/2)+2*beta+(-1)*IdA0*(1+(-1)*beta)**(1/2)*beta+alpha*(1+(-1)*beta)**(1/2)*beta+(-1)*beta**2+IdA0*(1+(-1)*beta)**(1/2)*beta**2+(-1)*IdA0*beta*gamma+IdA0*beta**2*gamma+gamma**2+alpha*(1+(-1)*beta)**(1/2)*gamma**2+(-1)*beta*gamma**2+(-1)*IaA0*(1+(1+(-1)*beta)**(1/2))*((-1)+beta)*((-1)+beta+gamma**2)+V0*(1+beta**2+(-1)*gamma**2+beta*((-2)+gamma**2))));
    #Idep=sym.exp(1)**(((-1)*t+t0)*gamma)*IdA0;





    t0_val=0
    #vtm=np.mean(vol_tr)/66.35
    #vrm=-1.0
    #d_in=st_point_dep1[r]
    #a_in=st_point_ada1[r]
    aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0,d_in)
    print('***** dep init ****')
    print(d_in)
    print('***** ada init ****')
    print(a_in)

    lam_x = sym.lambdify(t, aux, modules=['numpy'])
    x_vals = tim_aux / 15.58
    y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
    t_next=x_vals[y_vals>vtm][0]
    t_out.append(t_next)

    ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0, d_in).subs(t,t_next)
    dep = Idep.subs(gamma, gam).subs(t0, t0_val).subs(IdA0, d_in).subs(t,t_next)


    #with open('spt_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    spt = pickle.load(f)

    #with open('Iadap_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    Iadap_mat = pickle.load(f)

    #with open('Idep_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + 'gam_' + str(gam) + '_.pkl', "rb") as f:
    #    Idep_mat = pickle.load(f)

    #auxa = (IaA - a_in) ** 2
    #auxd = (IdA - d_in) ** 2
    #ada_nc=Iadap_mat[4,auxa.argmin(),auxd.argmin()]
    #dep_nc=Idep_mat[4, auxa.argmin(),auxd.argmin()]
    x_sector = np.logical_and(t0_val < x_vals, x_vals < t_next) * x_vals


    x_vals = x_sector[np.nonzero(x_sector)]
    y_vals = lam_x(x_vals)

    plt.figure()
    plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)


    #vrm=-0.7349960338647149
    #d_inc=dep_vec1_2[r,c]
    #a_inc=ada_vec1_2[r,c]
    for i in range(n_sp):
        t_init=t_next+2/15.58
        print('******dep ****'+str(i)+' *****spike****')
        print(dep + d_inc)
        print('***** ada ****'+str(i)+' *****spike****')
        print(ada + a_inc)


        aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        #print('ada comp')
        #print( ada+a_inc)
        #print('dep comp')
        #print(dep + d_inc)

#        print('ada')
#        print(ada_nc)
#        print('dep')
#        print(dep_nc)
        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / 15.58)>t_init]/ 15.58
        y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        t_init = t_next + 2 / 15.58

        if len(np.nonzero(y_vals > vtm)[0]) > 0:



            t_next=x_vals[y_vals>vtm][0]
            t_out.append(t_next)
            ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
            dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
#        auxa = (IaA - a_in) ** 2
#        auxd = (IdA - d_in) ** 2
#        ada_nc = Iadap_mat[4, auxa.argmin(), auxd.argmin()]
#        dep_nc = Idep_mat[4, auxa.argmin(), auxd.argmin()]
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            #plt.figure()
            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
        else:
            print('aoa')
            print(n_sp)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)

            i=n_sp+1

    plt.plot(tim, vol)
    return t_out


def plot_tr_inter(Vconvfact,vtm,vrm,a_in,d_in,a_inc,d_inc,bet,gam,n_sp,cell_num,alpha3,intervalli):
    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load

    tim_aux = np.linspace(0.0, intervalli[len(intervalli)-1][1], 10000)
    plt.figure()

    for ind_int in range(len(intervalli)):
        st_sign = intervalli[ind_int][0]
        len_sign = (intervalli[ind_int][1]-intervalli[ind_int][0])/15.58
        corr=int(alpha3/0.252587)*0.2


        #alpha3 = [0.252587 * 4, 0.252587 * 5]
        t_out = []
        tr=int(alpha3/0.252587)+10
        #cell_num = '95810005'


        t = sym.Symbol('t')
        alpha,beta,gamma,IaA0,IdA0,t0,V0= sym.symbols('alpha,beta,gamma,IaA0,IdA0,t0,V0')

        [V, Iadap, Idep] = load()


        if ind_int==0:
            t0_val = 0
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0,d_in)
        else:
            t0_val = t_next#(intervalli[ind_int][1]-intervalli[ind_int-1][1])/15.58
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, y_vals[len(y_vals) - 1]).subs(IaA0,a_in).subs(IdA0, d_in)

        print('***** dep init ****')
        print(d_in)
        print('***** ada init ****')
        print(a_in)

        lam_x = sym.lambdify(t, aux, modules=['numpy'])

        x_vals = tim_aux[(tim_aux / 15.58) > t0_val] / 15.58

        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        if len(x_vals[y_vals>vtm])>0:
            t_next=t0_val+min(x_vals[y_vals>vtm][0]-t0_val,len_sign)
        else:
            t_next=t0_val+len_sign
        t_out.append(t_next)

        ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0, d_in).subs(t,t_next)
        dep = Idep.subs(gamma, gam).subs(t0, t0_val).subs(IdA0, d_in).subs(t,t_next)


        x_sector = np.logical_and(t0_val < x_vals, x_vals < t_next) * x_vals


        x_vals = x_sector[np.nonzero(x_sector)]
        y_vals = lam_x(x_vals)

        #plt.figure()
        if ind_int==0:
            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
        else:
            plt.plot( x_vals * 15.58, y_vals * Vconvfact)

        stop=False
        i=0
        while(not(stop)):
            t_init = t_next + 2 / 15.58

            print('******dep ****'+str(i)+' *****spike****')
            print(dep + d_inc)
            print('***** ada ****'+str(i)+' *****spike****')
            print(ada + a_inc)
            if ind_int==0:
                aux1=0
            else:
                aux1=(st_sign/15.58)

            if t_init < aux1+len_sign:


                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
                lam_x = sym.lambdify(t, aux, modules=['numpy'])
                x_vals = tim_aux[(tim_aux / 15.58)>t_init]/ 15.58
                y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
                t_init = t_next + 2 / 15.58

                if len(x_vals[y_vals > vtm])>0:
                    t_next=t0_val+min(x_vals[y_vals>vtm][0],len_sign)
                else:
                    t_next =t0_val+ len_sign
                t_out.append(t_next)
                ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
                dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
                x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


                x_vals = x_sector[np.nonzero(x_sector)]
                y_vals = lam_x(x_vals)
            #plt.figure()

                if ind_int == 0:
                    plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
                else:
                    plt.plot(x_vals * 15.58, y_vals * Vconvfact)

                i=i+1
            else:
                #tim_aux = np.linspace(tim.min(), tim.max(), 10000)
                if (ind_int<len(intervalli)-1):
                    t_init=intervalli[ind_int][1]/15.58
                    t_next=intervalli[ind_int+1][0]/15.58
                    print('si riparte da................')
                    print(y_vals[len(y_vals)-1])
                    print('....................................')
                    aux = V.subs(alpha, 0.0).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, y_vals[len(y_vals)-1]).subs(IaA0,ada).subs(IdA0, dep)
                    ada = Iadap.subs(alpha, 0).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, y_vals[len(y_vals) - 1]).subs(IaA0, ada).subs(IdA0, dep).subs(t, t_next)
                    dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep).subs(t, t_next)

                    lam_x = sym.lambdify(t, aux, modules=['numpy'])
                    x_vals = tim_aux / 15.58
                    x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals

                    x_vals = x_sector[np.nonzero(x_sector)]

                    y_vals = lam_x(x_vals)
                # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
                    #t_next = (intervalli[ind_int+1][0]-st_sign)/15.58
                    t_out.append(t_next)


                    #plt.figure()
                    plt.plot( x_vals * 15.58, y_vals * Vconvfact)

                    print('aoa')
                    print(n_sp)
                    print(len(np.nonzero(y_vals > vtm)[0]))
                    print(vtm)
                i = i + 1
                stop = True

#        plt.plot(tim, vol)
    return t_out


def plot_tr_inter2(Vconvfact,vtm,vrm,a_in,d_in,a_inc,d_inc,bet,gam,n_sp,cell_num,alpha3,intervalli,tc):
    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load

    tim_aux = np.linspace(0.0, intervalli[len(intervalli)-1][1], 10000)
    plt.figure()
    t_out = []

    for ind_int in range(len(intervalli)):
        st_sign = intervalli[ind_int][0]
        len_sign = (intervalli[ind_int][1]-intervalli[ind_int][0])/15.58
        corr=int(alpha3/0.252587)*0.2


        #alpha3 = [0.252587 * 4, 0.252587 * 5]

        tr=int(alpha3/0.252587)+10
        #cell_num = '95810005'


        t = sym.Symbol('t')
        alpha,beta,gamma,IaA0,IdA0,t0,V0,ada0,dep0= sym.symbols('alpha,beta,gamma,IaA0,IdA0,t0,V0,ada0,dep0')

        [V, Iadap, Idep] = load()

        V2 = -1 + (-(V0) - 1) * sym.exp(1) ** (-300 * (t-t0))

        if ind_int==0:
            t0_val = 0
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0,d_in)
            #V2.subs(t,t0).subs(V0,-1) -Vconvfact + (-(V0 * Vconvfact) - Vconvfact) * sym.exp(1) ** (-300 * t)
            #aux=V2.subs(t0, t0_val).subs(V0, -1)
        else:
            t0_val = t_next#(intervalli[ind_int][1]-intervalli[ind_int-1][1])/15.58
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, y_vals[len(y_vals) - 1]).subs(IaA0,a_in).subs(IdA0, d_in)

        print('***** dep init ****')
        print(d_in)
        print('***** ada init ****')
        print(a_in)

        lam_x = sym.lambdify(t, aux, modules=['numpy'])

        x_vals = tim_aux[(tim_aux / 15.58) > t0_val] / 15.58

        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        if len(x_vals[y_vals>vtm])>0:
            t_next=t0_val+min(x_vals[y_vals>vtm][0]-t0_val,len_sign)
        else:
            t_next=t0_val+len_sign
        t_out.append(t_next)

        ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0, d_in).subs(t,t_next)
        dep = Idep.subs(gamma, gam).subs(t0, t0_val).subs(IdA0, d_in).subs(t,t_next)


        x_sector = np.logical_and(t0_val < x_vals, x_vals < t_next) * x_vals


        x_vals = x_sector[np.nonzero(x_sector)]
        y_vals = lam_x(x_vals)

        #plt.figure()
        if ind_int==0:
            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
        else:
            plt.plot( x_vals * 15.58, y_vals * Vconvfact)

        stop=False
        i=0
        while(not(stop)):
            t_init = t_next + 2 / 15.58

            print('******dep ****'+str(i)+' *****spike****')
            print(dep + d_inc)
            print('***** ada ****'+str(i)+' *****spike****')
            print(ada + a_inc)
            if ind_int==0:
                aux1=0
            else:
                aux1=(st_sign/15.58)

            if t_init < aux1+len_sign:


                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
                lam_x = sym.lambdify(t, aux, modules=['numpy'])
                x_vals = tim_aux[(tim_aux / 15.58)>t_init]/ 15.58
                y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
                t_init = t_next + 2 / 15.58

                if len(x_vals[y_vals > vtm])>0:
                    t_next=t0_val+min(x_vals[y_vals>vtm][0]-t0_val,len_sign)
                else:
                    t_next =t0_val+ len_sign
                t_out.append(t_next)
                ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
                dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
                x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


                x_vals = x_sector[np.nonzero(x_sector)]
                y_vals = lam_x(x_vals)
            #plt.figure()

                if ind_int == 0:
                    plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
                else:
                    plt.plot(x_vals * 15.58, y_vals * Vconvfact)

                i=i+1
            else:
                #tim_aux = np.linspace(tim.min(), tim.max(), 10000)
                if (ind_int<len(intervalli)-1):
                    t_init=intervalli[ind_int][1]/15.58
                    t_next=intervalli[ind_int+1][0]/15.58
                    print('si riparte da................')
                    print(y_vals[len(y_vals)-1])
                    print('....................................')
                    aux = V.subs(alpha, 0.0).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, y_vals[len(y_vals)-1]).subs(IaA0,ada).subs(IdA0, dep)
                    V2 = -1 + (V0 + 1) * sym.exp(1) ** (- (t - t0)/tc)
                    Iadap2=a_in+(ada0-a_in) * sym.exp(1) ** (- (t - t0)/tc)
                    Idep2=d_in+(dep0-d_in) * sym.exp(1) ** (- (t - t0)/tc)

                    aux = V2.subs(t0, t_init).subs(V0, y_vals[len(y_vals) - 1])
                    ada = Iadap2.subs(t0, t_init).subs(ada0,ada).subs(t, t_next)
                    dep = Idep2.subs(t0, t_init).subs(dep0, dep).subs(t, t_next)

                    #ada = Iadap.subs(alpha, 0).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, y_vals[len(y_vals) - 1]).subs(IaA0, ada).subs(IdA0, dep).subs(t, t_next)
                    #dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep).subs(t, t_next)

                    lam_x = sym.lambdify(t, aux, modules=['numpy'])
                    x_vals = tim_aux / 15.58
                    x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals

                    x_vals = x_sector[np.nonzero(x_sector)]
                    y_vals = np.zeros(len(x_vals))

                    for indic in range(len(x_vals)):
                        y_vals[indic]=aux.subs(t,x_vals[indic])


                # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
                    #t_next = (intervalli[ind_int+1][0]-st_sign)/15.58
                    t_out.append(t_next)


                    #plt.figure()
                    if (np.size(y_vals) == 1):
                        y_vals = np.ones(np.size(x_vals)) * y_vals

                    plt.plot( x_vals * 15.58, y_vals * Vconvfact)

                    print('aoa')
                    print(n_sp)
                    print(len(np.nonzero(y_vals > vtm)[0]))
                    print(vtm)
                i = i + 1
                stop = True

#        plt.plot(tim, vol)
    return t_out



def plot_tr_v3(Vconvfact,vtm,vrm,a_in,d_in,a_inc,d_inc,bet,delta1,n_sp,cell_num,alpha3,sc):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    st_sign = 531.1
    corr=alpha3/sc
    f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

    aux=f.readline()
    tim=[]
    vol=[]
    trovato = False
    for aux in f:
        tim.append(np.float64(aux.split('\t')[0]))
        vol.append(np.float64(aux.split('\t')[1]))
    tim = np.array(tim)
    vol = np.array(vol)

    tim_aux=np.linspace(tim.min(), tim.max(), 10000)

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
    aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t0_val).subs(V0, -1).subs(IaA0,a_in).subs(IdA0,d_in).subs(Psi, psi1)

    print('***** dep init ****')
    print(d_in)
    print('***** ada init ****')
    print(a_in)

    lam_x = sym.lambdify(t, aux, modules=['numpy'])
    x_vals = tim_aux / 15.58
    y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
    t_next=x_vals[y_vals>vtm][0]
    t_out.append(t_next)
################################################################################
    #ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0, d_in).subs(t,t_next)
    #dep = Idep.subs(gamma, gam).subs(t0, t0_val).subs(IdA0, d_in).subs(t,t_next)
#################################################################################

    ada=Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t0_val).subs(V0, -1).subs(IaA0,a_in).subs(IdA0, d_in).subs(Psi, psi1).subs(t, t_next)
    dep= Idep.subs(beta, bet).subs(IdA0, d_in).subs(t, t_next).subs(t0, t0_val)

    #with open('spt_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    spt = pickle.load(f)

    #with open('Iadap_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    Iadap_mat = pickle.load(f)

    #with open('Idep_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + 'gam_' + str(gam) + '_.pkl', "rb") as f:
    #    Idep_mat = pickle.load(f)

    #auxa = (IaA - a_in) ** 2
    #auxd = (IdA - d_in) ** 2
    #ada_nc=Iadap_mat[4,auxa.argmin(),auxd.argmin()]
    #dep_nc=Idep_mat[4, auxa.argmin(),auxd.argmin()]
    x_sector = np.logical_and(t0_val < x_vals, x_vals < t_next) * x_vals


    x_vals = x_sector[np.nonzero(x_sector)]
    y_vals = lam_x(x_vals)

    plt.figure()
    plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)


    #vrm=-0.7349960338647149
    #d_inc=dep_vec1_2[r,c]
    #a_inc=ada_vec1_2[r,c]
    for i in range(n_sp):
        t_init=t_next+2/15.58
        print('******dep ****'+str(i)+' *****spike****')
        print(dep + d_inc)
        print('***** ada ****'+str(i)+' *****spike****')
        print(ada + a_inc)


        #aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0, dep+d_inc).subs(Psi, psi1)

        #print('ada comp')
        #print( ada+a_inc)
        #print('dep comp')
        #print(dep + d_inc)

#        print('ada')
#        print(ada_nc)
#        print('dep')
#        print(dep_nc)
        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / 15.58)>t_init]/ 15.58
        y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        t_init = t_next + 2 / 15.58

        if len(np.nonzero(y_vals > vtm)[0]) > 0:



            t_next=x_vals[y_vals>vtm][0]
            t_out.append(t_next)
            ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(Psi, psi1).subs(t, t_next)
            dep = Idep.subs(beta, bet).subs(IdA0, dep+d_inc).subs(t, t_next).subs(t0, t_init)
            #ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
            #dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
#        auxa = (IaA - a_in) ** 2
#        auxd = (IdA - d_in) ** 2
#        ada_nc = Iadap_mat[4, auxa.argmin(), auxd.argmin()]
#        dep_nc = Idep_mat[4, auxa.argmin(), auxd.argmin()]
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            #plt.figure()
            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
        else:
            print('aoa')
            print(n_sp)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)

            i=n_sp+1

    plt.plot(tim, vol)
    return t_out


def plot_tr_v3_vect(Vconvfact,vtm,vrm,a_in,d_in,a_inc,d_inc,bet,delta1,alpha3,sc):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    n_sp=np.size(a_inc)
    st_sign = 531.1
    corr=alpha3/sc
    f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

    aux=f.readline()
    tim=[]
    vol=[]
    trovato = False
    for aux in f:
        tim.append(np.float64(aux.split('\t')[0]))
        vol.append(np.float64(aux.split('\t')[1]))
    tim = np.array(tim)
    vol = np.array(vol)

    tim_aux=np.linspace(tim.min(), tim.max(), 10000)

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
    aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t0_val).subs(V0, -1).subs(IaA0,a_in).subs(IdA0,d_in).subs(Psi, psi1)

    print('***** dep init ****')
    print(d_in)
    print('***** ada init ****')
    print(a_in)

    lam_x = sym.lambdify(t, aux, modules=['numpy'])
    x_vals = tim_aux / 15.58
    y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
    t_next=x_vals[y_vals>vtm][0]
    t_out.append(t_next)
################################################################################
    #ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t0_val).subs(V0, -1).subs(IaA0, a_in).subs(IdA0, d_in).subs(t,t_next)
    #dep = Idep.subs(gamma, gam).subs(t0, t0_val).subs(IdA0, d_in).subs(t,t_next)
#################################################################################

    ada=Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t0_val).subs(V0, -1).subs(IaA0,a_in).subs(IdA0, d_in).subs(Psi, psi1).subs(t, t_next)
    dep= Idep.subs(beta, bet).subs(IdA0, d_in).subs(t, t_next).subs(t0, t0_val)


    #with open('spt_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    spt = pickle.load(f)

    #with open('Iadap_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + '_v0_' + str(vrm) + 'bet_' + str(bet) + 'gam_' + str(gam) + '_.pkl',"rb") as f:
    #    Iadap_mat = pickle.load(f)

    #with open('Idep_surf_IaA_' + str(IaA_max) + '_' + str(IaA_min) + '_IdA_' + str(IdA_max) + '_' + str(IdA_min) + '_bin_' + str(bin_x_unt) + 'gam_' + str(gam) + '_.pkl', "rb") as f:
    #    Idep_mat = pickle.load(f)

    #auxa = (IaA - a_in) ** 2
    #auxd = (IdA - d_in) ** 2
    #ada_nc=Iadap_mat[4,auxa.argmin(),auxd.argmin()]
    #dep_nc=Idep_mat[4, auxa.argmin(),auxd.argmin()]
    x_sector = np.logical_and(t0_val < x_vals, x_vals < t_next) * x_vals


    x_vals = x_sector[np.nonzero(x_sector)]
    y_vals = lam_x(x_vals)

    plt.figure()
    plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)


    #vrm=-0.7349960338647149
    #d_inc=dep_vec1_2[r,c]
    #a_inc=ada_vec1_2[r,c]
    for i in range(n_sp):
        t_init=t_next+2/15.58
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('******dep ****'+str(i)+' *****spike****')
        print(dep )
        print('******dep inc ****' + str(i) + ' *****spike****')
        print(d_inc[i])
        print('******dep tot ****' + str(i) + ' *****spike****')
        print(dep + d_inc[i])

        #print('***** ada ****'+str(i)+' *****spike****')
        #print(ada)
        #print('***** ada inc ****' + str(i) + ' *****spike****')
        #print( a_inc[i])
        print('******ada tot ****' + str(i) + ' *****spike****')
        print(ada + a_inc[i])

        #aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, ada+a_inc[i]).subs(IdA0, dep+d_inc[i]).subs(Psi, psi1)

        #print('ada comp')
        #print( ada+a_inc)
        #print('dep comp')
        #print(dep + d_inc)

#        print('ada')
#        print(ada_nc)
#        print('dep')
#        print(dep_nc)
        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / 15.58)>t_init]/ 15.58
        y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        t_init = t_next + 2 / 15.58

        if len(np.nonzero(y_vals > vtm)[0]) > 0:



            t_next=x_vals[y_vals>vtm][0]
            t_out.append(t_next)
            print('***** t_next****')
            print(t_next)

            ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc[i]).subs(IdA0, dep+d_inc[i]).subs(Psi, psi1).subs(t, t_next)
            dep = Idep.subs(beta, bet).subs(IdA0, dep+d_inc[i]).subs(t, t_next).subs(t0, t_init)
            #ada = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0, t_init).subs(V0, vrm).subs(IaA0,ada+a_inc).subs(IdA0, dep+d_inc).subs(t, t_next)
            #dep = Idep.subs(gamma, gam).subs(t0, t_init).subs(IdA0, dep+d_inc).subs(t, t_next)
#        auxa = (IaA - a_in) ** 2
#        auxd = (IdA - d_in) ** 2
#        ada_nc = Iadap_mat[4, auxa.argmin(), auxd.argmin()]
#        dep_nc = Idep_mat[4, auxa.argmin(), auxd.argmin()]
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals


            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            #plt.figure()
            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
        else:
            print('aoa')
            print(n_sp)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)

            i=n_sp+1

    plt.plot(tim, vol)
    plt.title('trace (' + str(sc * alpha3) + 'nA)')
    return t_out


def plot_tr_v3_vect2(Vconvfact,vtm,vrm,a_inc,d_inc,bet,delta1,alpha3,sc):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    n_sp=np.size(a_inc)
    st_sign = 531.1
    corr=alpha3/sc
    f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

    aux=f.readline()
    tim=[]
    vol=[]
    trovato = False
    for aux in f:
        tim.append(np.float64(aux.split('\t')[0]))
        vol.append(np.float64(aux.split('\t')[1]))
    tim = np.array(tim)
    vol = np.array(vol)

    tim_aux=np.linspace(tim.min(), tim.max(), 10000)

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
    plt.figure()
    t_next=0
    for i in range(n_sp):
        t_init=t_next+2/15.58
        print('******dep inc ****' + str(i) + ' *****spike****')
        print(d_inc[i])

        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada inc ****' + str(i) + ' *****spike****')
        print( a_inc[i])

        #aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, a_inc[i]).subs(IdA0, d_inc[i]).subs(Psi, psi1)

        #print('ada comp')
        #print( ada+a_inc)
        #print('dep comp')
        #print(dep + d_inc)

#        print('ada')
#        print(ada_nc)
#        print('dep')
#        print(dep_nc)
        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / 15.58)>t_init]/ 15.58
        y_vals = lam_x(x_vals)
    # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)
        t_init = t_next + 2 / 15.58

        if len(np.nonzero(y_vals > vtm)[0]) > 0:



            t_next=x_vals[y_vals>vtm][0]
            t_out.append(t_next)

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

            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
        else:
            print('aoa')
            print(n_sp)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)
            t_next = 500/15.58
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals
            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)

            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)

            i=n_sp+1

    plt.plot(tim, vol)
    plt.title('trace ('+ str(sc*alpha3)+'nA)')
    return t_out



def plot_tr_v3_vect3_dep_ini(Vconvfact,vtm,vrm,a_inc,bet,delta1,alpha3,sc, tao_m, Idep_ini_vr, st_sign, dur_sign,Idep_ini):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    plotta=False
    n_sp=np.size(a_inc)
    corr=alpha3

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
    for i in range(n_sp+1):
        if i>0:
            t_init=t_next+2/tao_m
        else:
            t_init=t_next


        #aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if i==0:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0, Idep_ini).subs(Psi, psi1)
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
    return t_out,time,voltage

def plot_tr_v3_vect3_dep_ini2(Vconvfact,vtm,vrm,a_inc,bet,delta1,alpha3,sc, tao_m, Idep_ini_vr, st_sign, dur_sign,Idep_ini,Ith):

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
    for i in range(n_sp+1):
        if i>0:
            t_init=t_next+2/tao_m
        else:
            t_init=t_next


        #aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if i==0:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0, Idep_ini*(corr-Ith)*(corr>Ith)).subs(Psi, psi1)
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
    return t_out,time,voltage


def plot_tr_from_fit(Vconvfact, vtm, vrm, funzione,popt, bet, delta1,alpha3, sc):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    #n_sp=np.size(a_inc)
    st_sign = 531.1
    corr=alpha3/sc
    f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

    aux=f.readline()
    tim=[]
    vol=[]
    trovato = False
    for aux in f:
        tim.append(np.float64(aux.split('\t')[0]))
        vol.append(np.float64(aux.split('\t')[1]))
    tim = np.array(tim)
    vol = np.array(vol)

    tim_aux=np.linspace(tim.min(), tim.max(), 10000)

    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
    #cell_num = '95810005'

    t = sym.Symbol('t')
    delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load_v3()
    psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)

    t_init=0
    plt.figure()
    i=0
    while t_init*15.58<500:

    #while i < 20:
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada inc ****' + str(i) + ' *****spike****')
        print(funzione(t_init*15.58,*popt))
        print('t_init')
        print(t_init)
        print(t_init * 15.58)
        print('bet')
        print(bet)
        print('delta1')
        print(delta1)
        print('vrm')
        print(vrm)
        print('psi1')
        print(psi1)
        print(i)

        # aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if i==0:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,funzione(t_init * 15.58,*popt)).subs(IdA0, 0).subs(Psi, psi1)

        else:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, funzione(t_init*15.58,*popt[0])).subs(IdA0, 0).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / 15.58) > t_init] / 15.58
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
            t_out.append(t_next)
            print("t_next .....")
            print(t_next)
            print(t_next*15.58)



            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals

            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)

            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
            #st_sign=st_sign#+t_next*15.58
            print("st_sign .....")
            print(st_sign)
        else:
            print('aoa')
            print(i)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)
            t_next=500
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals
            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)

            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)

        i = i + 1
        print('t_next')
        print(t_next)
        t_init = t_next + 2 / 15.58

    plt.plot(tim, vol)
    plt.title('trace (' + str(alpha3/sc) + 'nA)')
    return t_out


def plot_tr_from_fit2(Vconvfact, vtm, vrm, exp_cum,popt,popt_dep, bet, delta1,alpha3, sc):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    #n_sp=np.size(a_inc)
    st_sign = 531.1
    corr=alpha3/sc
    try:
        f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

        aux=f.readline()
        tim=[]
        vol=[]
        trovato = False
        for aux in f:
            tim.append(np.float64(aux.split('\t')[0]))
            vol.append(np.float64(aux.split('\t')[1]))
        tim = np.array(tim)
        vol = np.array(vol)

        tim_aux=np.linspace(tim.min(), tim.max(), 10000)
    except:
        print('no trace')
        f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_1.0.soma.v.txt","r")
        aux = f.readline()
        tim = []
        for aux in f:
            tim.append(np.float64(aux.split('\t')[0]))
        tim = np.array(tim)

        tim_aux = np.linspace(tim.min(), tim.max(), 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
    #cell_num = '95810005'

    t = sym.Symbol('t')
    delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load_v3()
    psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)

    t_init=0
    plt.figure()
    i=0
    #while i<20:
    while t_init * 15.58 < 500:
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada  ****' + str(i) + ' *****spike****')
        print(exp_cum(t_init*15.58,popt[0],popt[1],popt[2]))
        print('***** dep  ****' + str(i) + ' *****spike****')
        print(exp_cum(t_init * 15.58, popt_dep[0], popt_dep[1],popt_dep[2]))

        popt_dep
        print('t_init')
        print(t_init)
        print(t_init * 15.58)


        # aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if i==0:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,exp_cum(t_init * 15.58,popt[0],popt[1],popt[2])).subs(IdA0, exp_cum(t_init * 15.58, popt_dep[0], popt_dep[1], popt_dep[2])).subs(Psi, psi1)

        else:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, exp_cum(t_init*15.58,popt[0],popt[1],popt[2])).subs(IdA0, exp_cum(t_init * 15.58, popt_dep[0], popt_dep[1], popt_dep[2])).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / 15.58) > t_init] / 15.58
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
            t_out.append(t_next)
            print("t_next .....")
            print(t_next)
            print(t_next*15.58)


            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals

            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)

            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
            #st_sign=st_sign#+t_next*15.58
            print("st_sign .....")
            print(st_sign)
        else:
            print('aoa')
            print(i)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)
            t_next = 500/15.58
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals
            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            try:
                plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
            except:
                plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact*np.ones(x_vals.size))

        i = i + 1
        print('t_next')
        print(t_next)
        t_init = t_next + 2 / 15.58
    try:
        plt.plot(tim, vol)
    except:
        print('no trace')
    plt.title('trace (' + str(alpha3/sc) + 'nA)')
    return t_out


def plot_tr_from_fit3(Vconvfact, vtm, vrm, exp_cum,popt, bet, delta1,alpha3, sc):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    #n_sp=np.size(a_inc)
    st_sign = 531.1
    corr=alpha3/sc
    try:
        f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

        aux=f.readline()
        tim=[]
        vol=[]
        trovato = False
        for aux in f:
            tim.append(np.float64(aux.split('\t')[0]))
            vol.append(np.float64(aux.split('\t')[1]))
        tim = np.array(tim)
        vol = np.array(vol)

        tim_aux=np.linspace(tim.min(), tim.max(), 10000)
    except:
        print('no trace')
        f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_1.0.soma.v.txt","r")
        aux = f.readline()
        tim = []
        for aux in f:
            tim.append(np.float64(aux.split('\t')[0]))
        tim = np.array(tim)

        tim_aux = np.linspace(tim.min(), tim.max(), 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
    #cell_num = '95810005'

    t = sym.Symbol('t')
    delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load_v3()
    psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)

    t_init=0
    plt.figure()
    i=0
    #while i<20:
    while t_init * 15.58 < 500:
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada  ****' + str(i) + ' *****spike****')
        print(exp_cum(t_init*15.58,popt[0],popt[1],popt[2]))

        print('t_init')
        print(t_init)
        print(t_init * 15.58)


        # aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if i==0:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,exp_cum(t_init * 15.58,popt[0],popt[1],popt[2])).subs(IdA0, 0).subs(Psi, psi1)

        else:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, exp_cum(t_init*15.58,popt[0],popt[1],popt[2])).subs(IdA0, 0).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / 15.58) > t_init] / 15.58
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
            t_out.append(t_next)
            print("t_next .....")
            print(t_next)
            print(t_next*15.58)


            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals

            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)

            plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
            #st_sign=st_sign#+t_next*15.58
            print("st_sign .....")
            print(st_sign)
        else:
            print('aoa')
            print(i)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)
            t_next = 500/15.58
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals
            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            try:
                plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact)
            except:
                plt.plot(st_sign + x_vals * 15.58, y_vals * Vconvfact*np.ones(x_vals.size))

        i = i + 1
        print('t_next')
        print(t_next)
        t_init = t_next + 2 / 15.58
    try:
        plt.plot(tim, vol)
    except:
        print('no trace')
    plt.title('trace (' + str(alpha3/sc) + 'nA)')
    return t_out


def plot_tr_from_fit4(Vconvfact, vtm, vrm, exp_cum,popt, bet, delta1,alpha3, sc,tao,Idep_ini_vr,st_sign,dur_sign):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    #n_sp=np.size(a_inc)

    corr=round(alpha3*sc)/1000
    try:
        f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

        aux=f.readline()
        tim=[]
        vol=[]
        trovato = False
        for aux in f:
            tim.append(np.float64(aux.split('\t')[0]))
            vol.append(np.float64(aux.split('\t')[1]))
        tim = np.array(tim)
        vol = np.array(vol)

        tim_aux=np.linspace(tim.min(), tim.max(), 100000)
    except:
        print('no trace')
        f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_1.0.soma.v.txt","r")
        aux = f.readline()
        tim = []
        for aux in f:
            tim.append(np.float64(aux.split('\t')[0]))
        tim = np.array(tim)

        tim_aux = np.linspace(tim.min(), tim.max(), 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
    #cell_num = '95810005'

    t = sym.Symbol('t')
    delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load_v3()
    psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)

    t_init=0
    plt.figure()
    i=0
    #while i<20:
    while t_init * tao < 400:
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada  ****' + str(i) + ' *****spike****')
        if i==0:
            print(exp_cum(t_init * tao, popt[0], popt[1], popt[2]))
        else:
            print(exp_cum(t_next*tao,popt[0],popt[1],popt[2]))

        print('t_init')
        print(t_init)
        print(t_init * tao)


        # aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if i==0:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0, 0).subs(Psi, psi1)

        else:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, exp_cum(t_next*tao,popt[0],popt[1],popt[2])).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / tao) > t_init] / tao
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
            if (t_next*tao<400):
                t_out.append(t_next*tao)
                print("t_next .....")
                print(t_next)
                print(t_next*tao)


            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals

            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)

            plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact)
            #st_sign=st_sign#+t_next*15.58
            print("st_sign .....")
            print(st_sign)
        else:
            print('aoa')
            print(i)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)
            t_next = 500/tao
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals
            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            try:
                plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact)
            except:
                plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact*np.ones(x_vals.size))

        i = i + 1
        print('t_next')
        print(t_next)
        t_init = t_next + 2 / tao
    try:
        plt.plot(tim, vol)
    except:
        print('no trace')
    plt.title('trace (' + str(alpha3*sc/1000) + 'nA)')
    return t_out

def plot_tr_from_fit5(Vconvfact, vtm, vrm, exp_cum,popt, bet, delta1,alpha3, sc,tao,Idep_ini_vr,st_sign,dur_sign,Neuron_tr,tempo_soglia):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3
    plotta=False
    #n_sp=np.size(a_inc)

    corr=round(alpha3*sc)/1000
    try:
        if Neuron_tr:
            f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

            aux=f.readline()
            tim=[]
            vol=[]
            trovato = False
            for aux in f:
                tim.append(np.float64(aux.split('\t')[0]))
                vol.append(np.float64(aux.split('\t')[1]))
            tim = np.array(tim)
            vol = np.array(vol)

            tim_aux=np.linspace(tim.min(), tim.max(), 100000)
        else:
            tr = int(round(alpha3*sc)/200) + 10
            cell_num = '95810010'
            abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
            abf.setSweep(tr)
            print('trrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr')
            print(tr)
            vol_base = 0#-62.86698717948718
            vol = abf.sweepY + vol_base
            tim=abf.sweepX * 1000
            tim_aux = np.linspace(tim.min(), tim.max(), 100000)
        #plt.plot(, abf.sweepY)

    except:
        print('no trace')
        f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_1.0.soma.v.txt","r")
        aux = f.readline()
        tim = []
        for aux in f:
            tim.append(np.float64(aux.split('\t')[0]))
        tim = np.array(tim)

        tim_aux = np.linspace(tim.min(), tim.max(), 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
    #cell_num = '95810005'

    t = sym.Symbol('t')
    delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load_v3()
    psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)

    t_init=0
    if plotta:
        plt.figure()
    i=0
    #while i<20:
    while t_init * tao < dur_sign:
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada  ****' + str(i) + ' *****spike****')
        if i==0:
            print(exp_cum(t_init * tao, popt[0], popt[1]))
        else:
            new_ada=exp_cum(t_next*tao,popt[0],popt[1])*(((t_next*tao)>tempo_soglia)*15+1)
            print(new_ada)
            print(exp_cum(t_next*tao,popt[0],popt[1]))

        print('t_init')
        print(t_init)
        print(t_init * tao)


        # aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if t_init * tao== dur_sign:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0,adaf).subs(IdA0, depf).subs(Psi, psi1)

        else:
            if i==0:
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0, 0).subs(Psi, psi1)

            else:
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, new_ada).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / tao) > t_init] / tao
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
            if i == 0:
                adaf = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, 0).subs(IdA0, 0).subs(Psi, psi1).subs(t, t_next)
                depf = Idep.subs(beta, bet).subs(IdA0, 0).subs(t0, t_init).subs(t, t_next)
            else:
                adaf = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0,vrm).subs(IaA0, new_ada).subs(IdA0, Idep_ini_vr).subs(Psi, psi1).subs(t, t_next)
                depf = Idep.subs(beta, bet).subs(IdA0, new_ada).subs(t0, t_init).subs(t, t_next)
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
            t_next = 500/tao
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals
            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            try:
                if plotta:
                    plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact)
            except:
                if plotta:
                    plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact*np.ones(x_vals.size))

        i = i + 1
        print('t_next')
        print(t_next)
        if t_next==dur_sign:
            adaf=Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0,new_ada ).subs(IdA0, Idep_ini_vr).subs(Psi, psi1).subs(t,dur_sign)
            depf=Idep.subs(beta, bet).subs(IdA0, new_ada).subs(t0, t_init).subs(t,dur_sign)
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
    return t_out

def plot_tr_from_fit6(Vconvfact, vtm, vrm, exp_cum,popt, bet, delta1,alpha3, sc,tao,Idep_ini_vr,st_sign,dur_sign,Neuron_tr,tempo_soglia):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3

    #n_sp=np.size(a_inc)

    corr=round(alpha3*sc)/1000
    try:
        if Neuron_tr:
            f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

            aux=f.readline()
            tim=[]
            vol=[]
            trovato = False
            for aux in f:
                tim.append(np.float64(aux.split('\t')[0]))
                vol.append(np.float64(aux.split('\t')[1]))
            tim = np.array(tim)
            vol = np.array(vol)

            tim_aux=np.linspace(tim.min(), tim.max(), 100000)
        else:
            tr = int(round(alpha3*sc)/200) + 10
            #cell_num = '95810010'
            #abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
            #abf.setSweep(tr)
            #print('trrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr')
            #print(tr)
            #vol_base = 0#-62.86698717948718
            #vol = abf.sweepY + vol_base
            #tim=abf.sweepX * 1000
            tim_aux = np.linspace(0, 1000, 10000)
        #plt.plot(, abf.sweepY)

    except:
        print('no trace')
        #f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_1.0.soma.v.txt","r")
        aux = f.readline()
        tim = []
        for aux in f:
            tim.append(np.float64(aux.split('\t')[0]))
        tim = np.array(tim)

        tim_aux = np.linspace(0, 1000, 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
    #cell_num = '95810005'

    t = sym.Symbol('t')
    delta, Psi, alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

    [V, Iadap, Idep] = load_v3()
    psi1 = ((-4) * bet + ((1 + delta1) ** 2)) ** (0.5)

    t_init=0
    plt.figure()
    i=0
    #while i<20:
    while t_init * tao < dur_sign:
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada  ****' + str(i) + ' *****spike****')
        if i==0:
            print(exp_cum(t_init * tao, popt[0], popt[1]))
        else:
            print(exp_cum(t_next*tao,popt[0],popt[1])*(((t_next*tao)>tempo_soglia)*15+1))
            print(exp_cum(t_next*tao,popt[0],popt[1]))

        print('t_init')
        print(t_init)
        print(t_init * tao)


        # aux = V.subs(alpha, alpha3).subs(beta, bet).subs(gamma, gam).subs(t0,t_init ).subs(V0, vrm).subs(IaA0, ada+a_inc).subs(IdA0,dep+d_inc)
        if i==0:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0, 0).subs(Psi, psi1)

        else:
            aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, exp_cum(t_next*tao,popt[0],popt[1])*(((t_next*tao)>tempo_soglia)*10+1)).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / tao) > t_init] / tao
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
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

            plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact)
            #st_sign=st_sign#+t_next*15.58
            print("st_sign .....")
            print(st_sign)
        else:
            print('aoa')
            print(i)
            print(len(np.nonzero(y_vals > vtm)[0]))
            print(vtm)
            t_next = 500/tao
            x_sector = np.logical_and(t_init < x_vals, x_vals < t_next) * x_vals
            x_vals = x_sector[np.nonzero(x_sector)]
            y_vals = lam_x(x_vals)
            try:
                plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact)
            except:
                plt.plot(st_sign + x_vals * tao, y_vals * Vconvfact*np.ones(x_vals.size))

        i = i + 1
        print('t_next')
        print(t_next)
        if t_next==dur_sign:
            t_init=t_next
        else:
            t_init = t_next + 2 / tao
    try:
        plt.plot(tim, vol)
    except:
        print('no trace')
    plt.title('trace (' + str(alpha3*sc/1000) + 'nA)')
    return t_out

def plot_tr_from_fit_neuron_dep_ini(Vconvfact, vtm, vrm, funzione,popt, bet, delta1,alpha3, sc,tao,Idep_ini_vr,st_sign,dur_sign,Neuron_tr,tempo_soglia,dep_ini):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3
    plotta=False
    #n_sp=np.size(a_inc)

    corr=round(alpha3*sc)/1000
    try:
        if Neuron_tr:
            f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

            aux=f.readline()
            tim=[]
            vol=[]
            trovato = False
            for aux in f:
                tim.append(np.float64(aux.split('\t')[0]))
                vol.append(np.float64(aux.split('\t')[1]))
            tim = np.array(tim)
            vol = np.array(vol)

            tim_aux=np.linspace(tim.min(), tim.max(), 100000)
        else:
            tr = int(round(alpha3*sc)/200) + 10
            cell_num = '95810010'
            abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
            abf.setSweep(tr)
            print('trrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr')
            print(tr)
            vol_base = 0#-62.86698717948718
            vol = abf.sweepY + vol_base
            tim=abf.sweepX * 1000
            tim_aux = np.linspace(tim.min(), tim.max(), 100000)
        #plt.plot(, abf.sweepY)

    except:
        print('no trace')
        #f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_1.0.soma.v.txt","r")
        #aux = f.readline()
        tim = []
        #for aux in f:
         #   tim.append(np.float64(aux.split('\t')[0]))
        #tim = np.array(tim)

        tim_aux = np.linspace(0, 1000, 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
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
    while t_init * tao < dur_sign:
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada  ****' + str(i) + ' *****spike****')
        if i==0:
            print(funzione(t_init * tao, *popt))
        else:
            new_ada=funzione(t_next*tao,*popt)*(((t_next*tao)>tempo_soglia)*15+1)
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
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0, dep_ini).subs(Psi, psi1)

            else:
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, new_ada).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / tao) > t_init] / tao
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
            if i == 0:
                adaf = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, 0).subs(IdA0, dep_ini).subs(Psi, psi1).subs(t, t_next)
                depf = Idep.subs(beta, bet).subs(IdA0, dep_ini).subs(t0, t_init).subs(t, t_next)
            else:
                adaf = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0,vrm).subs(IaA0, new_ada).subs(IdA0, Idep_ini_vr).subs(Psi, psi1).subs(t, t_next)
                depf = Idep.subs(beta, bet).subs(IdA0, new_ada).subs(t0, t_init).subs(t, t_next)
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
            depf=Idep.subs(beta, bet).subs(IdA0, new_ada).subs(t0, t_init).subs(t,dur_sign)
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

def plot_tr_from_fit_neuron_dep_ini2(Vconvfact, vtm, vrm, funzione,popt, bet, delta1,alpha3, sc,tao,Idep_ini_vr,st_sign,dur_sign,Neuron_tr,tempo_soglia,dep_ini,Ith):

    import pyabf
    import matplotlib.pyplot as plt
    import sympy as sym
    import numpy as np
    from load_eq import load_v3
    plotta=False
    #n_sp=np.size(a_inc)

    corr=round(alpha3*sc)/1000
    try:
        if Neuron_tr:
            f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

            aux=f.readline()
            tim=[]
            vol=[]
            trovato = False
            for aux in f:
                tim.append(np.float64(aux.split('\t')[0]))
                vol.append(np.float64(aux.split('\t')[1]))
            tim = np.array(tim)
            vol = np.array(vol)

            tim_aux=np.linspace(tim.min(), tim.max(), 100000)
        else:
            tr = int(round(alpha3*sc)/200) + 10
            cell_num = '95810010'
            abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
            abf.setSweep(tr)
            print('trrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr')
            print(tr)
            vol_base = 0#-62.86698717948718
            vol = abf.sweepY + vol_base
            tim=abf.sweepX * 1000
            tim_aux = np.linspace(tim.min(), tim.max(), 100000)
        #plt.plot(, abf.sweepY)

    except:
        print('no trace')
        #f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_1.0.soma.v.txt","r")
        #aux = f.readline()
        tim = []
        #for aux in f:
         #   tim.append(np.float64(aux.split('\t')[0]))
        #tim = np.array(tim)

        tim_aux = np.linspace(0, 1000, 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
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
    while t_init * tao < dur_sign:
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('***** ada  ****' + str(i) + ' *****spike****')
        if i==0:
            print(funzione(t_init * tao, *popt))
        else:
            new_ada=funzione(t_next*tao,*popt)*(((t_next*tao)>tempo_soglia)*15+1)
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
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0, dep_ini*(corr*1000-Ith)*(corr*1000>Ith)).subs(Psi, psi1)

            else:
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, new_ada).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / tao) > t_init] / tao
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
            if i == 0:
                adaf = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, 0).subs(IdA0, dep_ini*(corr*1000-Ith)*(corr*1000>Ith)).subs(Psi, psi1).subs(t, t_next)
                depf = Idep.subs(beta, bet).subs(IdA0, dep_ini*(corr*1000-Ith)*(corr*1000>Ith)).subs(t0, t_init).subs(t, t_next)
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





def plot_tr_from_fit_neuron_rette_dep_ini(Vconvfact, vtm, vrm, funzione,popt, bet, delta1,alpha3, sc,tao,Idep_ini_vr,st_sign,dur_sign,Neuron_tr,lin_func_inf,vinc_inf,lin_func_sup,vinc_sup,Idep_ini):

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
    if corr*1000<vinc_inf:
        dur_sign=min(dur_sign,lin_func_inf.subs(I,corr*1000))

    if corr*1000>vinc_sup:
        dur_sign=min(dur_sign,lin_func_sup.subs(I,corr*1000))
    try:
        if Neuron_tr:
            f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

            aux=f.readline()
            tim=[]
            vol=[]
            trovato = False
            for aux in f:
                tim.append(np.float64(aux.split('\t')[0]))
                vol.append(np.float64(aux.split('\t')[1]))
            tim = np.array(tim)
            vol = np.array(vol)

            tim_aux=np.linspace(tim.min(), tim.max(), 100000)
        else:
            tr = int(round(alpha3*sc)/200) + 10
            cell_num = '95810010'
            abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
            abf.setSweep(tr)
            print('trrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr')
            print(tr)
            vol_base = 0#-62.86698717948718
            vol = abf.sweepY + vol_base
            tim=abf.sweepX * 1000
            tim_aux = np.linspace(tim.min(), tim.max(), 100000)
        #plt.plot(, abf.sweepY)

    except:
        print('no trace')
        #f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_1.0.soma.v.txt","r")
        #aux = f.readline()
        tim = []
        #for aux in f:
         #   tim.append(np.float64(aux.split('\t')[0]))
        #tim = np.array(tim)

        tim_aux = np.linspace(0, 1000, 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
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
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, -1).subs(IaA0,0).subs(IdA0,Idep_ini).subs(Psi, psi1)

            else:
                aux = V.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, new_ada).subs(IdA0, Idep_ini_vr).subs(Psi, psi1)


        lam_x = sym.lambdify(t, aux, modules=['numpy'])
        x_vals = tim_aux[(tim_aux / tao) > t_init] / tao
        y_vals = lam_x(x_vals)
        # plt.plot(31.1 + x_vals * 15.58, y_vals * 66.35)


        if len(np.nonzero(y_vals > vtm)[0]) > 0:

            t_next = x_vals[y_vals > vtm][0]#+t_init
            if i == 0:
                adaf = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0, vrm).subs(IaA0, 0).subs(IdA0, Idep_ini).subs(Psi, psi1).subs(t, t_next)
                depf = Idep.subs(beta, bet).subs(IdA0, Idep_ini).subs(t0, t_init).subs(t, t_next)
            else:
                adaf = Iadap.subs(alpha, alpha3).subs(beta, bet).subs(delta, delta1).subs(t0, t_init).subs(V0,vrm).subs(IaA0, new_ada).subs(IdA0, Idep_ini_vr).subs(Psi, psi1).subs(t, t_next)
                depf = Idep.subs(beta, bet).subs(IdA0, new_ada).subs(t0, t_init).subs(t, t_next)
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
            depf=Idep.subs(beta, bet).subs(IdA0, new_ada).subs(t0, t_init).subs(t,dur_sign)
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

def plot_tr_from_fit_neuron_rette_dep_ini2(Vconvfact, vtm, vrm, funzione,popt, bet, delta1,alpha3, sc,tao,Idep_ini_vr,st_sign,dur_sign,Neuron_tr,lin_func_inf,vinc_inf,lin_func_sup,vinc_sup,Idep_ini,Ith):

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
    if corr*1000<vinc_inf:
        dur_sign=min(dur_sign,lin_func_inf.subs(I,corr*1000))

    if corr*1000>vinc_sup:
        dur_sign=min(dur_sign,lin_func_sup.subs(I,corr*1000))
    try:
        if Neuron_tr:
            f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_"+str(corr)[0:3]+".soma.v.txt", "r")
    # f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\CA1_pyr\\r_seed2_0\\step_" + str(corr)[0:3] + ".soma.v.txt","r")

            aux=f.readline()
            tim=[]
            vol=[]
            trovato = False
            for aux in f:
                tim.append(np.float64(aux.split('\t')[0]))
                vol.append(np.float64(aux.split('\t')[1]))
            tim = np.array(tim)
            vol = np.array(vol)

            tim_aux=np.linspace(tim.min(), tim.max(), 100000)
        else:
            tr = int(round(alpha3*sc)/200) + 10
            cell_num = '95810010'
            abf = pyabf.ABF('C:\\Users\\INA\\Downloads\\PC-cAC\\PC-cAC\\tutti\\' + cell_num + '.abf')
            abf.setSweep(tr)
            print('trrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr')
            print(tr)
            vol_base = 0#-62.86698717948718
            vol = abf.sweepY + vol_base
            tim=abf.sweepX * 1000
            tim_aux = np.linspace(tim.min(), tim.max(), 100000)
        #plt.plot(, abf.sweepY)

    except:
        print('no trace')
        #f = open("C:\\Users\\INA\\PycharmProjects\\neuron_model\\CA1_pyr2\\r_seed2_0\\step_1.0.soma.v.txt","r")
        #aux = f.readline()
        tim = []
        #for aux in f:
         #   tim.append(np.float64(aux.split('\t')[0]))
        #tim = np.array(tim)

        tim_aux = np.linspace(0, 1000, 10000)
    #alpha3 = [0.252587 * 4, 0.252587 * 5]
    t_out = []
    tr=int(alpha3/sc)+10
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
                depf = Idep.subs(beta, bet).subs(IdA0, new_ada).subs(t0, t_init).subs(t, t_next)
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
            depf=Idep.subs(beta, bet).subs(IdA0, new_ada).subs(t0, t_init).subs(t,dur_sign)
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

