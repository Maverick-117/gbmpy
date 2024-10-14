# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 22:15:01 2022

@author: jhvo9
"""
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 23:32:13 2021

Fix reprogramming level. Vary feedback types over the same graph.
now features a vary-over-mu-relaxation-rate loop
3 Nov: now features a vary-over-beta loop
22 Nov: now does variation over schedule type

Created on Wed Nov  3 20:41:24 2021

@author: jhvo9
"""
import matplotlib.pyplot as plt;
import numpy as np;
from os import makedirs, chdir, replace
from os.path import exists
import matplotlib;
from funciones import calc_EOT, calc_BED_Frac, dU_dt;
from numpy import sqrt
import matplotlib as mpl
from matplotlib import cm

class MplColorHelper:

  def __init__(self, cmap_name, start_val, stop_val):
    self.cmap_name = cmap_name
    self.cmap = plt.get_cmap(cmap_name)
    self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
    self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

  def get_rgb(self, val):
    return self.scalarMap.to_rgba(val)
def calc_eigs6(U,V):
    #assumes: using feedback type 6, with feedback on death, k2 model, l = 0, h1 = 1e5, h2 = 0, hd = h1/3, c = 5.169e-3
    eig1 = (1/(2*(1. + 100000.*V)**2*(123. + 4.e6*V)**2))*((-(1/(1. + 100000.*V)))*(-26.88877870433699 -7.126619396434031e6*V - 6.47097325564284e11*V**2 -2.3175998037183836e16*V**3 - 2.8436807407587528e20*V**4 -2151.1022963469572/(1. + 100000.*V) +(2.66198909172936e8*U)/(1. +100000.*V) - (8.429522819831166e8*V)/(1. + 100000.*V) +(7.055353202469684e13*U*V)/(1. +100000.*V) - (1.2178631508447506e14*V**2)/(1. + 100000.*V) +(6.406263523086407e18*U*V**2)/(1. +100000.*V) - (7.441343762417502e18*V**3)/(1. + 100000.*V) +(2.2944238056811978e23*U*V**3)/(1. +100000.*V) - (1.2809359896747807e23*V**4)/(1. + 100000.*V) +(2.815243933351163e27*U*V**4)/(1. +100000.*V) + (2.886335951870129e27*V**5)/(1. + 100000.*V) +(2.8436807407587454e31*V**6)/(1. + 100000.*V)) -sqrt((-(1/(1. + 100000.*V)**3))*4*(57840.51361686446 +3.799604471741991e10*V - 0.5*U*V +1.0674057901553622e16*V**2 +1.663339470839293e21*V**3 + 1.5543146926018233e26*V**4 +8.712766041868346e30*V**5 +2.686472848142641e35*V**6 + 3.014379014122504e39*V**7 -4.886936833640244e43*V**8 -1.6225927682921336e32*U*V**8 - 1.4798331884312886e48*V**9 -8.086520155362228e51*V**10) + (1/(1. + 100000.*V)**2)*(-26.88877870433699 -7.126619396434031e6*V -6.47097325564284e11*V**2 -2.3175998037183836e16*V**3 - 2.8436807407587528e20*V**4 -2151.1022963469572/(1. + 100000.*V) + (2.66198909172936e8*U)/(1. + 100000.*V) -(8.429522819831166e8*V)/(1. +100000.*V) + (7.055353202469684e13*U*V)/(1. +100000.*V) -(1.2178631508447506e14*V**2)/(1. +100000.*V) + (6.406263523086407e18*U*V**2)/(1. +100000.*V) -(7.441343762417502e18*V**3)/(1. +100000.*V) + (2.2944238056811978e23*U*V**3)/(1. +100000.*V) -(1.2809359896747807e23*V**4)/(1. +100000.*V) + (2.815243933351163e27*U*V**4)/(1. +100000.*V) +(2.886335951870129e27*V**5)/(1. +100000.*V) + (2.8436807407587454e31*V**6)/(1. +100000.*V))**2))
    eig2 = (1/(2*(1. + 100000.*V)**2*(123. + 4.e6*V)**2))*((-(1/(1. + 100000.*V)))*(-26.88877870433699 -7.126619396434031e6*V - 6.47097325564284e11*V**2 -2.3175998037183836e16*V**3 - 2.8436807407587528e20*V**4 -2151.1022963469572/(1. + 100000.*V) +(2.66198909172936e8*U)/(1. +100000.*V) - (8.429522819831166e8*V)/(1. + 100000.*V) +(7.055353202469684e13*U*V)/(1. +100000.*V) - (1.2178631508447506e14*V**2)/(1. + 100000.*V) +(6.406263523086407e18*U*V**2)/(1. +100000.*V) - (7.441343762417502e18*V**3)/(1. + 100000.*V) +(2.2944238056811978e23*U*V**3)/(1. +100000.*V) - (1.2809359896747807e23*V**4)/(1. + 100000.*V) +(2.815243933351163e27*U*V**4)/(1. +100000.*V) + (2.886335951870129e27*V**5)/(1. + 100000.*V) +(2.8436807407587454e31*V**6)/(1. + 100000.*V)) +sqrt((-(1/(1. + 100000.*V)**3))*4*(57840.51361686446 +3.799604471741991e10*V - 0.5*U*V +1.0674057901553622e16*V**2 +1.663339470839293e21*V**3 + 1.5543146926018233e26*V**4 +8.712766041868346e30*V**5 +2.686472848142641e35*V**6 + 3.014379014122504e39*V**7 -4.886936833640244e43*V**8 -1.6225927682921336e32*U*V**8 - 1.4798331884312886e48*V**9 -8.086520155362228e51*V**10) + (1/(1. + 100000.*V)**2)*(-26.88877870433699 -7.126619396434031e6*V -6.47097325564284e11*V**2 -2.3175998037183836e16*V**3 - 2.8436807407587528e20*V**4 -2151.1022963469572/(1. + 100000.*V) + (2.66198909172936e8*U)/(1. + 100000.*V) -(8.429522819831166e8*V)/(1. +100000.*V) + (7.055353202469684e13*U*V)/(1. +100000.*V) -(1.2178631508447506e14*V**2)/(1. +100000.*V) + (6.406263523086407e18*U*V**2)/(1. +100000.*V) -(7.441343762417502e18*V**3)/(1. +100000.*V) + (2.2944238056811978e23*U*V**3)/(1. +100000.*V) -(1.2809359896747807e23*V**4)/(1. +100000.*V) + (2.815243933351163e27*U*V**4)/(1. +100000.*V) +(2.886335951870129e27*V**5)/(1. +100000.*V) + (2.8436807407587454e31*V**6)/(1. +100000.*V))**2))
    return [eig1, eig2];

font = {'weight' : 'bold',
        'size'   : 100}

BASE_SIZE = 15;
SMALL_SIZE = 2.0 * BASE_SIZE;
MEDIUM_SIZE = 2.25 * BASE_SIZE;
LARGE_SIZE = 3 * BASE_SIZE;

mpl.rc('font', size=BASE_SIZE)          # controls default text sizes
mpl.rc('axes', titlesize=LARGE_SIZE)     # fontsize of the axes title
mpl.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
mpl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
mpl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
mpl.rc('legend', fontsize=BASE_SIZE)    # legend fontsize
mpl.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title

l_w = 10**(-7); # weak feedback on prob
l_s = 10**5; # strong feedback on prob
h1 = 10**5; # feedback on css div 
h2 = 10**5; # feedback on dcc div
l_vec =  [0,  0, l_w, l_s, l_w , l_s , 0 ,0 ,l_w,l_w,l_s,l_s];
h1_vec = [0, h1, 0  , 0  , h1  , h1  , h1,0 ,h1 ,0  , h1, 0 ];
h2_vec = [0, h2, 0  , 0  , h2  , h2  , 0 ,h2,0  ,h2 , 0 , h2];
h1a = 0;
#skip 10 20 40 45

#[1,6,7] had...the same stability behaviour/eigenvalue plots
saveQ = False;
cases = ["2 Gy", "2.67 Gy" , "3.4 Gy", "5 Gy"]; 
comp_str = "Gy";
Doses = [2.0,2.0];
Fracs = [calc_BED_Frac(.17,.02,Doses[0]),calc_BED_Frac(.17,.02,Doses[1])]#list(map(lambda d:calc_BED_Frac(a,b,d,40*(1+(3.4)/8.5)),Doses));#;[25,20]
EOTs  = [calc_EOT(Fracs[0],100), calc_EOT(Fracs[1],100)];
EOT_diff = EOTs[0] - EOTs[1];
EOT_str = ",EOT-shifted"
match_EOT =True;
if not match_EOT:
    EOT_diff = 0.;
    EOT_str = "";
pwr = 3
C = [5.196*10**(-pwr), 5.196*10**(-pwr)];
comp_dir = str(C[0])+"_reprog";
fix_str= str(Doses[-1])+"Gy"
comp_list = C;
unit = ' Gy';
T_stop = 1000;
#lll_load = 6;
step = 1;

# ss
load_temp = [3];#np.arange(12).tolist();#[1,3,5,6,7,10,11];#np.arange(12).tolist();
#deathFdbkQ = True;
death_r2_mult = 1.1;
log_yscaleQ = False;
if log_yscaleQ:
    log_str = " (log)";
    log_dir = "\\log_yscale\\";
else:
    log_str = "";
    log_dir = "\\linear_yscale\\"

for deathFdbkQ in [False]:    
    if deathFdbkQ:
        hd = 10**0;#32520.32520325203 * 10 ** (np.log10(3));#1e5; 
        deathFdbk_str = "_w_death_fdbk";
    else:
        hd = 0.0;#1e5; 
        deathFdbk_str = "_w_no_death_fdbk";
    kimReprogQ = False;
    kimDynamQ = False;
    if kimReprogQ:
        sub_drty = "kim_reprog";
    else:
        sub_drty = "\\corrected_reprog";
    
    if kimDynamQ:
        dyn_str = "\\kim_dynam";
    else:
        dyn_str = "\\new_dynam";
    
    deathVal_dir = '\\death_val_of_kim';
    color = ['r','k','g','b','m'];
    
    DT = 3.9;
    r1 = np.log(2)/DT; # growth rate of CSC
    r2 = np.log(2)/DT; # growth rate of DCC
    p = .505; # self renewal probability 
    l_w = 10**(-7); # weak feedback on prob
    l_s = 10**5; # strong feedback on prob
    h1 = 10**5; # feedback on css div 
    h2 = 10**(5); # feedback on dcc div
    d = r2;
    dsc = 0;
    
    compBaselineQ = False;
    showReprogQ = True;
    showNonReprogQ = True;
    schedChangeQ = False;
    use_muQ = True;
    use_insetsQ = False;
    if use_muQ:
        #, 1e3,np.infty]
        # sig_vec = [10];#list(map(lambda p: 1 * 10 ** (p), np.arange(0,3).tolist()));
        # rho_vec = [0];#[1/5];#sorted(list(map(lambda p: 2 ** (p), np.arange(0,-5,-1).tolist())) + [1/5]);
        sig_vec = [1];#sorted(list(map(lambda p:  10 ** (p), np.arange(-1,3).tolist()))+[0]);
        rho_vec = [0];#0.2
    # rename rho to rho
        xi1_vec = [1];#[0.01];#sorted(list(map(lambda p:  10 ** (p), np.arange(-1,3).tolist()))+[0]);
        xi2_vec = [0.1];#[1];#sorted(list(map(lambda p:  10 ** (p), np.arange(-1,3).tolist()))+[0]);#np.arange(0.01,0.1,0.01);#[.1];#sorted(list(map(lambda p:  10 ** (p), np.arange(0,3).tolist()))+[0])
        xi3_vec = [1e9]; 
        s_xi1_vec = ['\\' + str(el) for el in xi1_vec];
        s_xi2_vec = ['\\' + str(elb) for elb in xi2_vec];
        s_xi3_vec = ['\\' + str(els) for els in xi3_vec];
        s_sig_vec = ['\\' + str(el) for el in sig_vec];
        s_rho_vec = ['\\' + str(elb) + '\\' for elb in rho_vec];
        case = 'final\\4_June';#'hU\\6_Jan';#'reversal\\21_Oct';#'mathCheck\\2_Feb'#'reversionAttempt4\\14_Nov';#'muNN-l2b\\2_Nov';#'baseline\\31_Oct';
    else:
        sig_vec=[0];
        s_sig_vec = ['\\0\\'];
        case = 'reprog_check\\5_May';#'full_fdbk_gains\\13_Oct'
    
    if schedChangeQ:
        Nsched_vec = [1];
    else:
        Nsched_vec = [0];
    # Nsched_vec = [0,1];
    Nsched_vec_str = ['',' (delayed)']
    
    goldenLaptopQ = True;
    
    if goldenLaptopQ:
        prefix = "C:\\Users\\jhvo9\\Documents\\";
    else:
        prefix = "G:\\My Drive";
    
    base_dir = prefix+"\\a PhD Projects\\GBM Modeling\\python scripts\\data\\k2_model\\"
    base_dir1 = base_dir+"conventional_BED\\"+case+deathVal_dir+"\\"+sub_drty + dyn_str;#"conventional_BED\\srv_test\\16_Oct\\death_val_of_kim\\corrected_reprog";
    base_dir2 = base_dir+"fixed_BED\\40\\03_Aug\\death_val_of_kim\\"+sub_drty + dyn_str;
    
    save_dir = "G:\My Drive\a PhD Projects\GBM Modeling\python scripts\plots\k2_model\long term steady state sims\12000 days\death feedback"
    
    
    data2 = np.ones((len(Nsched_vec),len(xi1_vec),len(xi2_vec),len(xi3_vec),len(rho_vec),len(sig_vec),12)).tolist();
    data1 = np.ones((len(Nsched_vec),len(xi1_vec),len(xi2_vec),len(xi3_vec),len(rho_vec),len(sig_vec),12)).tolist();
    dataN = np.ones((len(Nsched_vec),len(xi1_vec),len(xi2_vec),len(xi3_vec),len(rho_vec),len(sig_vec),12)).tolist();
    for lll_load in load_temp:
        for nn,eln in enumerate(Nsched_vec):
            s_fr_0 = '_'+str(Fracs[0])+'_Fracs_';
            s_fr_1 = '_'+str(Fracs[1])+'_Fracs_';#'_'+str(Fracs[1])+'_Fracs_';
    
            for ss, els in enumerate(s_xi3_vec):
                for tt, elt in enumerate(s_xi1_vec):
                    for uu, elu in enumerate(s_xi2_vec):
                        for ww, elb in enumerate(s_rho_vec):
                            for vv, el in enumerate(s_sig_vec):
                                s_mu_str = elt + elu + els + elb  + el +"\\Schedule"+str(eln) + "\\";
                                plot_str2 = "TU_"+str(C[0])+'_reprog'+str(Doses[0])+'_'+comp_str+s_fr_0+str(T_stop)+'days_'+deathFdbk_str+"_w_reprog_"+str(lll_load);
                                if showNonReprogQ:
                                    plot_str1 = "TU_"+str(C[1])+'_reprog'+str(Doses[1])+'_'+comp_str+s_fr_1+str(T_stop)+'days_'+deathFdbk_str+"_w_no_reprog_"+str(lll_load);
                                else:
                                    plot_str1 = "TU_"+str(C[1])+'_reprog'+str(Doses[1])+'_'+comp_str+s_fr_1+str(T_stop)+'days_'+deathFdbk_str+"_w_reprog_"+str(lll_load);
                                data2[nn][tt][uu][ss][ww][vv][lll_load] = np.loadtxt(base_dir1 + s_mu_str + plot_str2 + ".txt");
                                data1[nn][tt][uu][ss][ww][vv][lll_load] = np.loadtxt(base_dir1 + s_mu_str + plot_str1 + ".txt");
                                print('dir2',base_dir1 + s_mu_str + plot_str2 + ".txt")
                                print('dir1',base_dir1 + s_mu_str + plot_str1 + ".txt")
                                plotn_str2 = "TU_none_"+str(C[1])+'_reprog'+str(Doses[0])+'_'+comp_str+s_fr_0+str(T_stop)+'days_'+deathFdbk_str+"_w_reprog_"+str(lll_load); 
                                print('dir_none',base_dir1 + s_mu_str + plotn_str2 + ".txt")
                                dataN[nn][tt][uu][ss][ww][vv][lll_load] = np.loadtxt(base_dir1 + s_mu_str + plotn_str2 + ".txt");
    
    #T_stop_inset = 200;
    
    #T_stop = 15000;
    #T_stop = 2000;
    T_stop_inset = 200;
    #T_stop = 1000;
    total_cell_num = 4/3*np.pi*(1 ** 3)*1e9;
    s_sched0 = str((Doses[0],Fracs[0]));
    s_sched1 = str((Doses[1],Fracs[1]));
    s_reprg0 = f"{float(C[0]):.3}";
    s_reprg1 = f"{float(C[1]):.3}";
    s_reprg_title = ";c:";
    
    
    if use_muQ:
        width, height = 20, 20;# 7: 0.15, 0.15
        x1, x2 = 0.25, 0.75#6: 0.1, 0.6   7: 0.3, 0.8
        y1, y2 = 0.65, 0.2 #6: 0.65, 0.2    7: 0.57, 0.05
        ticksize = 11;
    
        s_sched_title = ";sched:"
        s_title_tweakable = s_sched_title+s_sched0;
        s_lgd_0 = s_sched0;#str(csr[0])+unit+EOT_str;
        s_lgd_1 = s_sched1+";no reprog;"+EOT_str;#str(csr[1])+unit;
        COL2 = MplColorHelper('viridis', 0, len(load_temp)-1)#len(Nsched_vec)-1)
        COL1 = MplColorHelper('cool', 0, len(load_temp)-1)
        COLN = MplColorHelper('autumn', 0, len(load_temp)-1)
        for lll_load in load_temp:
            s_p = ""; s_ru = ""; s_rv = ""; s_d = "";
            h1 = h1_vec[lll_load];
            h2 = h2_vec[lll_load];
            l = l_vec[lll_load];
            figCellNo, axCellNo = plt.subplots(ncols=2, nrows=2, figsize=(width, height),
                                constrained_layout=True)
            #figCellNo, axCellNo[0,0] = plt.subplots(figsize=(8,8));
            #figCellNo, axCellNo[0,1] = plt.subplots(figsize=(8,8));
            #s_title = "(l:"+f"{float(l_vec[lll_load]):.3}"+";h1:"+f"{float(h1_vec[lll_load]):.3}"+";h2:"+f"{float(h2_vec[lll_load]):.3}"+";hd:"+f"{hd:.3}"+s_title_tweakable+")";
            if l_vec[lll_load] > 0:
                s_p = "p, ";
            if h1_vec[lll_load] > 0:
                s_ru = "r_u, ";
            if h2_vec[lll_load] > 0:
                s_rv = "r_v, ";
            if deathFdbkQ:
                s_d = "d";
            if l_vec[lll_load]+h1_vec[lll_load]+h2_vec[lll_load]+hd > 0:
                s_title = "Feedback on " +s_p + s_ru + s_rv + s_d;
            else:
                s_title = "Baseline Model";
            # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
            # s_lgd_1 = "w/ constant reprog";#"Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
    
            figCellNo.suptitle(s_title)
            for nn,eln in enumerate(Nsched_vec):
                # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
                # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
                for tt, elt in enumerate(xi1_vec):
                    xi1 = elt;
                    s_xi1 = '('+str(elt)+',';
                    for uu, elu in enumerate(xi2_vec):
                        xi2 = elu;
                        s_xi2 = str(elu)+',';
                        for ss, els in enumerate(xi3_vec):
                            xi3 = els;
                            s_xi3 = str(els)+',';
                            for ww, elb in enumerate(rho_vec):
                                beta = elb;
                                s_beta = str(elb) + ',';
                                for vv, el in enumerate(sig_vec):
                                    s_sig = '(σ='+str(el) +')';
                                    s_title_tweakable = s_xi1+s_xi2+s_beta+s_sig;
                                    st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                    if showReprogQ:
                                        axCellNo[0,0].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step] * total_cell_num,'--',color=COL2.get_rgb(nn),label='CSC '+s_lgd_0+'_'+str(lll_load))
                                        axCellNo[0,1].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step] * total_cell_num,'-',color=COL2.get_rgb(nn),label='DCC '+s_lgd_0+s_sig+'_'+str(lll_load))
                                    if showNonReprogQ:
                                        st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                        axCellNo[0,0].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step] * total_cell_num,'--',color=COL1.get_rgb(nn),label='CSC '+s_lgd_1+s_sig+'_'+str(lll_load))
                                        axCellNo[0,1].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step] * total_cell_num,'-',color=COL1.get_rgb(nn),label='DCC '+s_lgd_1+s_sig+'_'+str(lll_load))
                            
            st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
            axCellNo[0,0].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step] * total_cell_num,'g--',label='CSC, no treatment_'+str(lll_load))
            axCellNo[0,1].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step] * total_cell_num,'g-',label='DCC, no treatment_'+str(lll_load))
            axCellNo[0,0];#.legend(fancybox=True, framealpha=0.5);
            axCellNo[0,0].set_ylabel("CSC Cell Number"+log_str)#,fontsize=2*ticksize)
            axCellNo[0,0].set_xlabel("time (days)")#,fontsize=2*ticksize)
            #axCellNo[0,0].set_title(s_title)#,fontsize=2.5*ticksize)
            #axCellNo[0,0].tick_params(axis='both', which='major',labelsize=1.5*ticksize)
            axCellNo[0,1];#.legend(fancybox=True, framealpha=0.5);
            axCellNo[0,1].set_ylabel("DCC Cell Number"+log_str)#,fontsize=2*ticksize)
            axCellNo[0,1].set_xlabel("time (days)")#,fontsize=2*ticksize)
            #axCellNo[0,1].set_title(s_title)#,fontsize=2.5*ticksize)
            #axCellNo[0,1].tick_params(axis='both', which='major',labelsize=1.5*ticksize)
            if log_yscaleQ:
                axCellNo[0,0].set_yscale('log')
                axCellNo[0,1].set_yscale('log')
            if saveQ:
                figCellNo.savefig(plot_drty+"pop"+plot_save_suffix+".png",dpi=300)
                    #replace("pop"+plot_save_suffix+".png",plot_drty+"pop"+plot_save_suffix+".png")
            if use_insetsQ:
                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                axCellNo00 = figCellNo.add_axes([x1, y1, width, height])
                if showReprogQ:
                    axCellNo00.plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step] * total_cell_num,'--',color=COL2.get_rgb(nn),label='CSC '+s_lgd_0+'_'+str(lll_load))
                if showNonReprogQ:
                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                    axCellNo00.plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step]+EOT_diff,data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step] * total_cell_num,'--',color=COL1.get_rgb(nn),label='CSC '+s_lgd_1+s_sig+'_'+str(lll_load))
                st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max()+1;
                axCellNo00.plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step] * total_cell_num,'g--',label='CSC, no treatment_'+str(lll_load))
                axCellNo00.set_yscale('log')
                axCellNo00.tick_params(axis='both', which='major',labelsize=ticksize)
        
        
        
                axCellNo01 = figCellNo.add_axes([x2, y1, width, height]);
                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                if showReprogQ:
                    axCellNo01.plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step] * total_cell_num,'--',color=COL2.get_rgb(nn),label='DCC '+s_lgd_0+'_'+str(lll_load))
        
                if showNonReprogQ:
                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                    axCellNo01.plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step]+EOT_diff,data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step] * total_cell_num,'-.',color=COL1.get_rgb(nn),label='DCC '+s_lgd_1+s_sig+'_'+str(lll_load))
                st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max()+1;
                axCellNo01.plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop+1],dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop+1] * total_cell_num,'g--',label='DCC, no treatment_'+str(lll_load))
                axCellNo01.set_yscale('log')
                axCellNo01.tick_params(axis='both', which='major',labelsize=ticksize)
            
                # figCellNo
        
                axCellNo10 = figCellNo.add_axes([x1, y2, width, height]);
                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                if showReprogQ:
                     axCellNo10.plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],(data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]+data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]) * total_cell_num,'--',color=COL2.get_rgb(nn),label='DCC '+s_lgd_0+'_'+str(lll_load))
        
                if showNonReprogQ:
                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                    axCellNo10.plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step]+EOT_diff,(data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]+data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]) * total_cell_num,'-.',color=COL1.get_rgb(nn),label='DCC '+s_lgd_1+s_sig+'_'+str(lll_load))
                st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max()+1;
                axCellNo10.plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],(dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]+dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]) * total_cell_num,'g--',label='DCC, no treatment_'+str(lll_load))
                axCellNo10.set_yscale('log')
                axCellNo10.tick_params(axis='both', which='major')
            # if saveQ:
            #     figCellNo.savefig("total_pop"+plot_save_suffix+".png",dpi=300)
            #     # replace(base_dirGD+"total_pop"+plot_save_suffix+".png",plot_drty+"total_pop"+plot_save_suffix+".png")
                     #figCellNo, axTo = plt.subplots(figsize=(8,8));
            for nn,eln in enumerate(Nsched_vec):
                # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
                # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
                for tt, elt in enumerate(xi1_vec):
                    xi1 = elt;
                    s_xi1 = '\\'+str(elt)+'\\';
                    for uu, elu in enumerate(xi2_vec):
                        xi2 = elu;
                        s_xi2 = '\\'+str(elu)+'\\';
                        for ww, elb in enumerate(rho_vec):
                            beta = elb;
                            for vv, el in enumerate(sig_vec):
                                s_sig = '(σ='+ str(el) +')';
                                
                                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                
                                if log_yscaleQ:
                                    axCellNo[1,0].set_yscale('log')
                                if showReprogQ:
                                    axCellNo[1,0].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],(data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]+data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]) * total_cell_num,'--',color=COL2.get_rgb(nn),label='total '+s_lgd_0+s_sig+'_'+str(lll_load))
                                axCellNo[1,0].set_ylabel("Total Cell Number"+log_str)#,fontsize=2*ticksize)
                                axCellNo[1,0].set_xlabel("time (days)")#,fontsize=2*ticksize)
                                #axCellNo[1,0].set_title(s_title)#,fontsize=2.5*ticksize)
                                #axCellNo[1,0].tick_params(axis='both', which='major',labelsize=1.5*ticksize)
                                if showNonReprogQ:
                                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                    axCellNo[1,0].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,(data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]+data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]) * total_cell_num,'--',color=COL1.get_rgb(nn),label='total '+s_lgd_1+s_sig+'_'+str(lll_load))
                            st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                            axCellNo[1,0].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],(dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]+dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]) * total_cell_num,color[2]+'--',label='total, no treatment_'+str(lll_load))
            axCellNo[1,0];#.legend(fancybox=True, framealpha=0.5);
        # TODO: maybe add the other ones too?
            #figC, axC = plt.subplots(figsize=(8,8));
            for nn,eln in enumerate(Nsched_vec):
                # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
                # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
                for tt, elt in enumerate(xi1_vec):
                    xi1 = elt;
                    s_xi1 = '\\'+str(elt)+'\\';
                    for uu, elu in enumerate(xi2_vec):
                        xi2 = elu;
                        s_xi2 = '\\'+str(elu)+'\\';
                        for ww, elb in enumerate(rho_vec):
                            beta = elb;
                            for vv, el in enumerate(sig_vec):
                                s_sig = '(σ='+ str(el) +')';
                                
                                csc_p = [];
                                csc_p += [data2[nn][tt][uu][ss][ww][vv][lll_load][1]/(data2[nn][tt][uu][ss][ww][vv][lll_load][1] + data2[nn][tt][uu][ss][ww][vv][lll_load][2])];
                                csc_p += [data1[nn][tt][uu][ss][ww][vv][lll_load][1]/(data1[nn][tt][uu][ss][ww][vv][lll_load][1] + data1[nn][tt][uu][ss][ww][vv][lll_load][2])];
                                
                                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();

                                
                                if log_yscaleQ:
                                    axCellNo[1,1].set_yscale('log')
                                if showReprogQ:
                                    axCellNo[1,1].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],csc_p[0][:st2_stop:step],'-',color=COL2.get_rgb(nn),label='CSC frac'+s_lgd_0+s_sig+'_'+str(lll_load))
                                
                                csc_p_inset = [];
                                csc_p_inset += [data2[nn][tt][uu][ss][ww][vv][lll_load][1]/(data2[nn][tt][uu][ss][ww][vv][lll_load][1] + data2[nn][tt][uu][ss][ww][vv][lll_load][2])];
                                csc_p_inset += [data1[nn][tt][uu][ss][ww][vv][lll_load][1]/(data1[nn][tt][uu][ss][ww][vv][lll_load][1] + data1[nn][tt][uu][ss][ww][vv][lll_load][2])];
                                # for i in range(len(csc_p)):
                                #     for j in range(100000,len(csc_p[0])):
                                #         csc_p[i][j] = 0;
                                #         csc_p_inset[i][j] = 0;                                
                                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                                
                                if use_insetsQ:
                                    axCellNo11 = figCellNo.add_axes([x2, y2, width, height]);
                                    axCellNo11.plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],csc_p_inset[0][:st2_stop:step],'--',color=COL2.get_rgb(nn),label='DCC '+s_lgd_0+'_'+str(lll_load))
    
                                if showNonReprogQ:
                                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                    axCellNo[1,1].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,csc_p[1][:st1_stop:step],'-',color=COL1.get_rgb(nn),label='CSC frac'+s_lgd_1+s_sig+'_'+str(lll_load))
                                    if use_insetsQ:
                                        st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                                        axCellNo11.plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,csc_p[1][:st1_stop:step],'-',color=COL1.get_rgb(nn),label='CSC frac'+s_lgd_1+s_sig+'_'+str(lll_load))
                                
                # replace(base_dirGD+"csc"+plot_save_suffix+".png",plot_drty+"csc"+plot_save_suffix+".png")
            csc_p += [dataN[nn][tt][uu][ss][ww][vv][lll_load][1]/(dataN[nn][tt][uu][ss][ww][vv][lll_load][1] + dataN[nn][tt][uu][ss][ww][vv][lll_load][2])];
            csc_p_inset += [dataN[nn][tt][uu][ss][ww][vv][lll_load][1]/(dataN[nn][tt][uu][ss][ww][vv][lll_load][1] + dataN[nn][tt][uu][ss][ww][vv][lll_load][2])];
            st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
            # for i in range(len(csc_p)):
            #     for j in range(100000,len(csc_p[0])):
            #         csc_p[i][j] = 0;
            #         csc_p_inset[i][j] = 0;            
            axCellNo[1,1].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],csc_p[2][:st0_stop:step],color[2]+'-',label='CSC frac, no treatment')
            axCellNo[1,1];#.legend(fancybox=True, framealpha=0.5);
            axCellNo[1,1].set_ylabel("CSC frac")#,fontsize=2*ticksize)
            axCellNo[1,1].set_xlabel("time (days)")#),fontsize=2*ticksize)
            #axCellNo[1,1].set_title(s_title)#,fontsize=2.5*ticksize)
            #axCellNo[1,1].tick_params(axis='both', which='major',labelsize=1.5*ticksize)

                    
            if use_insetsQ:
                st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max()+1;
                axCellNo11.plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],csc_p_inset[2][:st0_stop:step],color[2]+'-',label='CSC frac, no treatment')
                axCellNo11.tick_params(axis='both', which='major',labelsize=ticksize)
            plt.show()
            # if saveQ:
            #     figC.savefig("csc"+plot_save_suffix+".png",dpi=300)
            figRest, axRest = plt.subplots(ncols=2, nrows=2, figsize=(2*width, height),
                                    constrained_layout=True)
            figRest.suptitle(s_title)
            
            #figTs, axTs = plt.subplots(figsize=(8,8));
            for nn,eln in enumerate(Nsched_vec):
                # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
                # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
                for tt, elt in enumerate(xi1_vec):
                    xi1 = elt;
                    s_xi1 = '\\'+str(elt)+'\\';
                    for uu, elu in enumerate(xi2_vec):
                        xi2 = elu;
                        s_xi2 = '\\'+str(elu)+'\\';
                        for ww, elb in enumerate(s_rho_vec):
                            beta = elb;
                            for vv, el in enumerate(sig_vec):
                                s_sig = '(σ='+ str(el) +')';
                                
                                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                
                                # if log_yscaleQ:
                                #     axRest[0,0].set_yscale('log')
                                if showReprogQ:                            
                                    axRest[0,0].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],data2[nn][tt][uu][ss][ww][vv][lll_load][3,:st2_stop:step],'--',color=COL2.get_rgb(nn),label='mu frac'+s_lgd_0+s_sig+'_'+str(lll_load))
                                if showNonReprogQ:
                                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                    axRest[0,0].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,data1[nn][tt][uu][ss][ww][vv][lll_load][3,:st1_stop:step],'--',color=COL1.get_rgb(nn),label='mu frac'+s_lgd_1+s_sig+'_'+str(lll_load))
    
            st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
            axRest[0,0].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],dataN[nn][tt][uu][ss][ww][vv][lll_load][3,:st0_stop:step],'g--',label='mu frac, no treatment_'+str(lll_load))
            axRest[0,0].legend(fancybox=True, framealpha=0.5);
                    
            axRest[0,0].set_ylabel("mu frac"+log_str)
            axRest[0,0].set_xlabel("time (days)")
            #axRest[0,0].set_title(s_title)
            if showNonReprogQ:
                st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                axRest[0,0].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,data1[nn][tt][uu][ss][ww][vv][lll_load][3,:st1_stop:step],'--',color=COL1.get_rgb(nn),label='mu frac'+s_lgd_1+s_sig+'_'+str(lll_load))
            
                # if saveQ:
                #     figTs.savefig(plot_drty+"pop"+plot_save_suffix+".png",dpi=300)
            #figTsr, axTsr = plt.subplots(figsize=(8,8));
            for nn,eln in enumerate(Nsched_vec):
                # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
                # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
                for tt, elt in enumerate(xi1_vec):
                    xi1 = elt;
                    s_xi1 = '\\'+str(elt)+'\\';
                    for uu, elu in enumerate(xi2_vec):
                        xi2 = elu;
                        s_xi2 = '\\'+str(elu)+'\\';
                        for ww, elb in enumerate(rho_vec):
                            beta = elb;
                            for vv, el in enumerate(sig_vec):
                                s_sig = '(σ='+ str(el) +')';
                                
                                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                if showReprogQ:                            
                                    axRest[0,1].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],beta * data2[nn][tt][uu][ss][ww][vv][lll_load][3,:st2_stop:step] * data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step] * total_cell_num,'--',color=COL2.get_rgb(nn),label='mu rate'+s_lgd_0+s_sig+'_'+str(lll_load))
                                axRest[0,1].set_ylabel("mu rate"+log_str)
                                axRest[0,1].set_xlabel("time (days)")
                                #axRest[0,1].set_title(s_title)
                                if showNonReprogQ:
                                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                    axRest[0,1].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,beta * data1[nn][tt][uu][ss][ww][vv][lll_load][3,:st1_stop:step] * data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step] * total_cell_num,'--',color=COL1.get_rgb(nn),label='mu rate'+s_lgd_1+s_sig+'_'+str(lll_load))
            
            st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
            axRest[0,1].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],beta * dataN[nn][tt][uu][ss][ww][vv][lll_load][3,:st0_stop:step] * dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step] * total_cell_num,'g--',label='mu rate, no treatment_'+str(lll_load))
            axRest[0,1].legend(fancybox=True, framealpha=0.5);
            # if saveQ:
            #     figTsr.savefig(plot_drty+"pop"+plot_save_suffix+".png",dpi=300)
        
        #xi3*U[0]/(1+xi3*U[0])
            #figTsrU, axTsrU = plt.subplots(figsize=(8,8));
            for nn,eln in enumerate(Nsched_vec):
                # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
                # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
                for tt, elt in enumerate(xi1_vec):
                    xi1 = elt;
                    s_xi1 = '\\'+str(elt)+'\\';
                    for uu, elu in enumerate(xi2_vec):
                        xi2 = elu;
                        s_xi2 = '\\'+str(elu)+'\\';
                        for ww, elb in enumerate(rho_vec):
                            beta = elb;
                            for vv, el in enumerate(sig_vec):
                                s_sig = '(σ='+ str(el) +')';
                                
                                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                if showReprogQ:                            
                                    axRest[1,0].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],xi3* data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]/(1+xi3* data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]),'--',color=COL2.get_rgb(nn),label='feedforward rate coefficient'+s_lgd_0+s_sig+'_'+str(lll_load))
                                axRest[1,0].set_ylabel("feedforward rate coefficient"+log_str)
                                axRest[1,0].set_xlabel("time (days)")
                                #axRest[1,0].set_title(s_title)
                                if showNonReprogQ:
                                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                    axRest[1,0].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,xi3 * data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]/(1+ xi3*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]),'--',color=COL1.get_rgb(nn),label='feedforward rate coefficient'+s_lgd_1+s_sig+'_'+str(lll_load))
            
            st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
            axRest[1,0].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],xi3 * dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]/(1+xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]),'g--',label='feedforward rate coefficient, no treatment_'+str(lll_load))
            axRest[1,0].legend(fancybox=True, framealpha=0.5);
    
            #figTsrVd, axTsrVd = plt.subplots(figsize=(8,8));
            for nn,eln in enumerate(Nsched_vec):
                # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
                # s_lgd_1 = "Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';
                for tt, elt in enumerate(xi1_vec):
                    xi1 = elt;
                    s_xi1 = '\\'+str(elt)+'\\';
                    for uu, elu in enumerate(xi2_vec):
                        xi2 = elu;
                        s_xi2 = '\\'+str(elu)+'\\';
                        for ww, elb in enumerate(rho_vec):
                            beta = elb;
                            for vv, el in enumerate(sig_vec):
                                s_sig = '(σ='+ str(el) +')';
                                
                                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                #
                                if showReprogQ:
                                    axRest[1,1].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],r2/(1+h2*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]) * xi3*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]/(1+xi3*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]) - d-(death_r2_mult*r2-d)*hd*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]/(1+hd*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]),'--',color=COL2.get_rgb(nn),label='effective division rate coefficient'+s_lgd_0+s_sig+'_'+str(lll_load))
                                axRest[1,1].set_ylabel("effective division rate coefficient"+log_str)
                                axRest[1,1].set_xlabel("time (days)")
                                #axRest[1,1].set_title(s_title)
                                if showNonReprogQ:
                                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                                    axRest[1,1].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,r2/(1+h2*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]) * xi3*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]/(1+xi3*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]) - d-(death_r2_mult*r2-d)*hd*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]/(1+hd*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]),'--',color=COL1.get_rgb(nn),label='effective division rate coefficient'+s_lgd_1+s_sig+'_'+str(lll_load))
            
            st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
            axRest[1,1].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],r2/(1+h2*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]) * xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]/(1+xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]) - d-(death_r2_mult*r2-d)*hd*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]/(1+hd*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]),'g--',label='effective division rate coefficient, no treatment_'+str(lll_load))
            axRest[1,1].legend(fancybox=True, framealpha=0.5);
            
            if use_insetsQ:
                axRest00_inset = figRest.add_axes([x1,y1, width, height])
                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                if showReprogQ:
                    axRest00_inset.plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],data2[nn][tt][uu][ss][ww][vv][lll_load][3,:st2_stop:step],'--',color=COL2.get_rgb(nn),label='mu frac'+s_lgd_0+s_sig+'_'+str(lll_load))
                if showNonReprogQ:
                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                    axRest00_inset.plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,data1[nn][tt][uu][ss][ww][vv][lll_load][3,:st1_stop:step],'--',color=COL1.get_rgb(nn),label='mu frac'+s_lgd_1+s_sig+'_'+str(lll_load))
                st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max()+1;
                axRest00_inset.plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],dataN[nn][tt][uu][ss][ww][vv][lll_load][3,:st0_stop:step],'g--',label='mu frac, no treatment_'+str(lll_load))
                axRest01_inset = figRest.add_axes([x2,y1, width, height])
                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                if showReprogQ:
                    axRest01_inset.plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],beta * data2[nn][tt][uu][ss][ww][vv][lll_load][3,:st2_stop:step] * data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step] * total_cell_num,'--',color=COL2.get_rgb(nn),label='mu rate'+s_lgd_0+s_sig+'_'+str(lll_load))
                if showNonReprogQ:
                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                    axRest01_inset.plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,beta * data1[nn][tt][uu][ss][ww][vv][lll_load][3,:st1_stop:step] * data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step] * total_cell_num,'--',color=COL1.get_rgb(nn),label='mu rate'+s_lgd_1+s_sig+'_'+str(lll_load))
                st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max()+1;
                axRest01_inset.plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],beta * dataN[nn][tt][uu][ss][ww][vv][lll_load][3,:st0_stop:step] * dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step] * total_cell_num,'g--',label='mu rate, no treatment_'+str(lll_load))
                
                axRest10_inset = figRest.add_axes([x1,y2, width, height])
                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                if showReprogQ:
                    axRest10_inset.plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],xi3* data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]/(1+xi3* data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]),'--',color=COL2.get_rgb(nn),label='feedforward rate coefficient'+s_lgd_0+s_sig+'_'+str(lll_load))
                if showNonReprogQ:
                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                    axRest10_inset.plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,xi3 * data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]/(1+ xi3*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]),'--',color=COL1.get_rgb(nn),label='feedforward rate coefficient'+s_lgd_1+s_sig+'_'+str(lll_load))
                st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max()+1;
                axRest10_inset.plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],xi3 * dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]/(1+xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]),'g--',label='feedforward rate coefficient, no treatment_'+str(lll_load))
                
                axRest11_inset = figRest.add_axes([x2,y2, width, height])
                st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                if showReprogQ:
                    axRest11_inset.plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],r2/(1+h2*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]) * xi3*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]/(1+xi3*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]) - d-(death_r2_mult*r2-d)*hd*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]/(1+hd*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]),'--',color=COL2.get_rgb(nn),label='effective division rate coefficient'+s_lgd_0+s_sig+'_'+str(lll_load))
                if showNonReprogQ:
                    st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max();
                    axRest11_inset.plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step]+EOT_diff,r2/(1+h2*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]) * xi3*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]/(1+xi3*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]) - d-(death_r2_mult*r2-d)*hd*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]/(1+hd*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]),'--',color=COL1.get_rgb(nn),label='effective division rate coefficient'+s_lgd_1+s_sig+'_'+str(lll_load))
                st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop_inset).max()+1;
                axRest11_inset.plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],r2/(1+h2*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]) * xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]/(1+xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]) - d-(death_r2_mult*r2-d)*hd*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]/(1+hd*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]),'g--',label='effective division rate coefficient, no treatment_'+str(lll_load))
            plt.show()
            
            figFdbk, axFdbk = plt.subplots(ncols=2, nrows=1, figsize=(width, width/2),
                                    constrained_layout=True)
            figDftn, axDftn = plt.subplots(ncols=2, nrows=1, figsize=(width, width/2),
                                    constrained_layout=True)
            figNetDiv, axNetDiv = plt.subplots(ncols=1, nrows=1, figsize=(width/2, width/2),
                                    constrained_layout=True)
            
            st2_stop = np.flatnonzero(data2[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
            if showReprogQ:
                axFdbk[0].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],r2/(1+h2*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]) * xi3*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]/(1+xi3*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]),'--',color=COL2.get_rgb(nn),label='DCC division rate coefficient_'+s_lgd_0+s_sig+'_'+str(lll_load))
                axFdbk[0].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],d+(death_r2_mult*r2-d)*hd*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]/(1+hd*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]),'-',color=COL2.get_rgb(nn),label='DCC death rate coefficient_'+s_lgd_0+s_sig+'_'+str(lll_load))
                axFdbk[1].plot(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],r1/(1+h1*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]/(1+h1a*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step])),'-',color=COL2.get_rgb(nn),label='CSC total division rate coefficient_'+s_lgd_0+s_sig+'_'+str(lll_load))
        
                axDftn[0].semilogy(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],(2*p/(1+l*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step])-1)*(r1/(1+h1*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]/(1+h1a*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step])) - dsc)*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]* total_cell_num,'--',color=COL2.get_rgb(nn),label='CSC self-renewal rate_'+s_lgd_0+s_sig+'_'+str(lll_load))
                axDftn[1].semilogy(data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],2*(1-p/(1+l*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]))*(r1/(1+h1*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]/(1+h1a*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step])) - dsc)*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]* total_cell_num,'-',color=COL2.get_rgb(nn),label='CSC differentiation rate_'+s_lgd_0+s_sig+'_'+str(lll_load))
            
                axNetDiv.semilogy(
                    data2[nn][tt][uu][ss][ww][vv][lll_load][0,:st2_stop:step],(2*(1-p/(1+l*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]))*(r1/(1+h1*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]/(1+h1a*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step])))*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step] + 
                                      (r2/(1+h2*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]) * xi3*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]/(1+xi3*data2[nn][tt][uu][ss][ww][vv][lll_load][1,:st2_stop:step]) - 
                                       d-(death_r2_mult*r2-d)*hd*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step]/(1+hd*data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step])) * data2[nn][tt][uu][ss][ww][vv][lll_load][2,:st2_stop:step])* total_cell_num,'-',color=COL2.get_rgb(nn),label='DCC overall division_'+s_lgd_0+s_sig+'_'+str(lll_load))
            if showNonReprogQ:                                                                                
                st1_stop = np.flatnonzero(data1[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
                axFdbk[0].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step],r2/(1+h2*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]),'--',color=COL1.get_rgb(nn),label='DCC division rate coefficient_'+s_lgd_1+s_sig+'_'+str(lll_load))
                axFdbk[0].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step],d+(death_r2_mult*r2-d)*hd*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]/(1+hd*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]),'-',color=COL1.get_rgb(nn),label='DCC death rate coefficient_'+s_lgd_1+s_sig+'_'+str(lll_load))
                axFdbk[1].plot(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step],r1/(1+h1*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]/(1+h1a*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step])),'-',color=COL1.get_rgb(nn),label='CSC total division rate coefficient_'+s_lgd_1+s_sig+'_'+str(lll_load))
        
                axDftn[0].semilogy(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step],(2*p/(1+l*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step])-1)*(r1/(1+h1*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]/(1+h1a*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step])))*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]* total_cell_num,'--',color=COL1.get_rgb(nn),label='CSC self-renewal rate_'+s_lgd_1+s_sig+'_'+str(lll_load))
                axDftn[1].semilogy(data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step],2*(1-p/(1+l*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]))*(r1/(1+h1*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]/(1+h1a*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step])))*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]* total_cell_num,'-',color=COL1.get_rgb(nn),label='CSC differentiation rate_'+s_lgd_1+s_sig+'_'+str(lll_load))
        
                axNetDiv.semilogy(
                    data1[nn][tt][uu][ss][ww][vv][lll_load][0,:st1_stop:step],(2*(1-p/(1+l*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]))*(r1/(1+h1*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]/(1+h1a*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step])))*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step] + 
                                      (r2/(1+h2*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]) * xi3*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]/(1+xi3*data1[nn][tt][uu][ss][ww][vv][lll_load][1,:st1_stop:step]) - 
                                       d-(death_r2_mult*r2-d)*hd*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]/(1+hd*data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step])) * data1[nn][tt][uu][ss][ww][vv][lll_load][2,:st1_stop:step]
                                      )* total_cell_num,'-',color=COL1.get_rgb(nn),label='DCC overall division_'+s_lgd_1+s_sig+'_'+str(lll_load))
        
    
            st0_stop = np.flatnonzero(dataN[nn][tt][uu][ss][ww][vv][lll_load][0] < T_stop).max();
            axFdbk[0].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],r2/(1+h2*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]) * xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]/(1+xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]),'g--',label='DCC division rate coefficient, no treatment_'+str(lll_load))
            axFdbk[0].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],d+(death_r2_mult*r2-d)*hd*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]/(1+hd*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]),'g-',label='DCC death rate coefficient, no treatment_'+str(lll_load))
            axFdbk[1].plot(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],r1/(1+h1*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]/(1+h1a*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step])),'g-',label='CSC total division rate coefficient, no treatment_'+str(lll_load))
    
            axDftn[0].semilogy(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],(2*p/(1+l*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step])-1)*(r1/(1+h1*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]/(1+h1a*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]))) * dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]* total_cell_num,'g--',label='CSC self-renewal rate, no treatment_'+str(lll_load))
            axDftn[1].semilogy(dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],2*(1-p/(1+l*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]))*(r1/(1+h1*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]/(1+h1a*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]))) * dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]* total_cell_num,'g-',label='CSC differentiation rate, no treatment_'+str(lll_load))
    
            axNetDiv.semilogy(
                dataN[nn][tt][uu][ss][ww][vv][lll_load][0,:st0_stop:step],(2*(1-p/(1+l*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]))*(r1/(1+h1*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]/(1+h1a*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step])))*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step] + 
                                  (r2/(1+h2*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]) * xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]/(1+xi3*dataN[nn][tt][uu][ss][ww][vv][lll_load][1,:st0_stop:step]) - 
                                   d-(death_r2_mult*r2-d)*hd*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]/(1+hd*dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step])) * dataN[nn][tt][uu][ss][ww][vv][lll_load][2,:st0_stop:step]
                                  )* total_cell_num,'g-',label='DCC overall division, no treatment_'+str(lll_load))
    
            #axFdbk[0,0].axhline(y = r2, color = 'r', linestyle = '-')
            axFdbk[0].legend(fancybox=True, framealpha=0.5);
            axFdbk[1].legend(fancybox=True, framealpha=0.5);
            axDftn[0].legend(fancybox=True, framealpha=0.5);
            axDftn[1].legend(fancybox=True, framealpha=0.5);
            axFdbk[0].set_ylabel("DCC rate coefficients"+log_str)
            axFdbk[0].set_xlabel("time (days)")
            axFdbk[1].set_xlabel("time (days)")
            axFdbk[1].set_ylabel("CSC total division rate coefficients"+log_str)
            axFdbk[1].set_xlabel("time (days)")
            axDftn[0].set_ylabel("CSC self-renewal rate"+log_str)
            axDftn[1].set_xlabel("time (days)")
            axDftn[1].set_ylabel("CSC differentiation rate"+log_str)
            axDftn[1].set_xlabel("time (days)")
            axDftn[0].set_ylabel('CSC renewal rate'+log_str)
            axDftn[0].set_xlabel('time (days)')
            axNetDiv.set_ylabel("Net DCC division rate")
            axNetDiv.set_xlabel("time (days)")
            axNetDiv.legend(fancybox=True, framealpha=0.5);
            figFdbk.suptitle(s_title)
            figDftn.suptitle(s_title)
            axNetDiv.set_title(s_title)
            plt.show()
            plt.show()
            plt.show()
            #axFdbk[0,0].axhline(y = r2, color = 'r', linestyle = '-')
            # figDCC, axDCC = plt.subplots(figsize=(8,8))
            # axDCC.
            
    
    if compBaselineQ:
        s_fr_0 = '_'+str(Fracs[0])+'_Fracs_';
        s_fr_1 = '_'+str(Fracs[1])+'_Fracs_';#'_'+str(Fracs[1])+'_Fracs_';
        base_dir = prefix+"\\a PhD Projects\\GBM Modeling\\python scripts\\data\\k2_model\\"
        caseC = "baseline_beta\\5_Nov"
        base_dirC = base_dir+"conventional_BED\\"+caseC+"\\death_val_of_kim\\"+sub_drty;
        data2C = np.ones((len(rho_vec),len(sig_vec),12)).tolist();
        data1C = np.ones((len(rho_vec),len(sig_vec),12)).tolist();
        dataNC = np.ones((len(rho_vec),len(sig_vec),12)).tolist();
        T_stop = 500;
        rho_vec = [0];
        sig_vec = [0];
        for lll_load in load_temp:
            for ww, elb in enumerate(rho_vec):
                for vv, el in enumerate(sig_vec):
                    plot_str2 = "TU_"+str(C[0])+'_reprog'+str(Doses[0])+'_'+comp_str+s_fr_0+str(T_stop)+'days_'+deathFdbk_str+"_"+str(lll_load);
                    plot_str1 = "TU_"+str(C[1])+'_reprog'+str(Doses[1])+'_'+comp_str+s_fr_1+str(T_stop)+'days_'+deathFdbk_str+"_"+str(lll_load);
                    data2C[ww][vv][lll_load] = np.loadtxt(base_dirC + "\\0\\0" + "\\" + plot_str2 + ".txt");
                    data1C[ww][vv][lll_load] = np.loadtxt(base_dirC + "\\0\\0" + "\\" + plot_str1 + ".txt");
                    plotn_str2 = "TU_none_"+str(C[0])+'_reprog'+str(Doses[0])+'_'+comp_str+s_fr_0+str(T_stop)+'days_'+deathFdbk_str+"_"+str(lll_load); 
        figVec = [figCellNo, figCellNo, figTs, figCellNo, figC, figTsr];
        axVec = [axCellNo[0,0], axCellNo[0,1], axTs, axTo, axC, axTsr];
        sig_vec = [1];
        for ww, elb in enumerate(rho_vec):
            for vv, el in enumerate(sig_vec):
                stC_stop = np.flatnonzero(data2C[ww][vv][lll_load][0] < T_stop).max();
                s_sig = '('+ str(el) +')';
                s_null = ', RT only no mu'
                axVec[0].plot(data2C[ww][vv][lll_load][0,:stC_stop:step],data2C[ww][vv][lll_load][1,:stC_stop:step] * total_cell_num,'c--',label='CSC'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
                axVec[1].plot(data2C[ww][vv][lll_load][0,:stC_stop:step],data2C[ww][vv][lll_load][2,:stC_stop:step] * total_cell_num,'c-',label='DCC'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
                axVec[2].plot(data2C[ww][vv][lll_load][0,:stC_stop:step],data2C[ww][vv][lll_load][3,:stC_stop:step],'c:',label='mu'+s_null+s_sig+'_'+str(lll_load))
                axVec[3].plot(data2C[ww][vv][lll_load][0,:stC_stop:step],(data2C[ww][vv][lll_load][1,:stC_stop:step]+data2C[ww][vv][lll_load][2,:stC_stop:step]) * total_cell_num,'c--',label='total'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
                csc_p_0 = data2C[ww][vv][lll_load][1]/(data2C[ww][vv][lll_load][1] + data2C[ww][vv][lll_load][2])
                axVec[4].plot(data2C[ww][vv][lll_load][0,:stC_stop:step],csc_p_0[:stC_stop:step],'c-',label='CSC'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
                axVec[5].plot(data2C[ww][vv][lll_load][0,:stC_stop:step],beta * data2C[ww][vv][lll_load][3,:stC_stop:step] * data2C[ww][vv][lll_load][2,:stC_stop:step] * total_cell_num,'c-',label='mu rate'+s_null+s_lgd_0+s_sig+'_'+str(lll_load))
        #[axVec[i].legend() for i in [0,1,2,3,4,5]]
        [axVec[i].set_xlim(0,T_stop+100) for i in [0,1,2,3,4,5]]
    # only relevant for feedback regime 6, model 0 or 1.
    # eigs2 = np.array(list(map(lambda a,b: calc_eigs6(a,b), data2[vv][lll_load][1,:],data2[vv][lll_load][2,:])))
    # eigs1 = np.array(list(map(lambda a,b: calc_eigs6(a,b), data1[vv][lll_load][1,:],data1[vv][lll_load][2,:])))
    #eigs0 = np.array(list(map(lambda a,b: calc_eigs6(a,b), dataN[lll_load][1,:],dataN[lll_load][2,:])))
    
    # for lll_load in load_temp:
    #     figE, axE = plt.subplots(figsize=(8,8));
    #     for vv, el in enumerate(sig_vec):
    #         s_sig = '('+ str(el[0])+','+str(el[1])+')';
        
    #         s_title = "(l:"+f"{float(l_vec[lll_load]):.3}"+";h1:"+f"{float(h1_vec[lll_load]):.3}"+";h2:"+f"{float(h2_vec[lll_load]):.3}"+";hd:"+f"{hd:.3}"+s_title_tweakable+")";
            
            
    #         st2_stop = np.flatnonzero(data2[vv][lll_load][0] < T_stop).max();
    #         st1_stop = np.flatnonzero(data1[vv][lll_load][0] < T_stop).max();
    #         # st0_stop = np.flatnonzero(dataN[lll_load][0] < T_stop).max();
        
    #         axE.plot(data2[vv][lll_load][0,:st2_stop:50],eigs2[:st2_stop:50,0],'_',color=COL2.get_rgb(vv),label='eig0 '+s_lgd_0+s_sig+'_'+str(lll_load))
    #         axE.plot(data2[vv][lll_load][0,:st2_stop:50],eigs2[:st2_stop:50,1],'|',color=COL2.get_rgb(vv),label='eig1 '+s_lgd_0+s_sig+'_'+str(lll_load))
    #         axE.plot(data1[vv][lll_load][0,:st1_stop:50],eigs1[:st1_stop:50,0],color[1]+'_',label='eig0 '+s_lgd_1)
    #         axE.plot(data1[vv][lll_load][0,:st1_stop:50],eigs1[:st1_stop:50,1],color[1]+'|',label='eig1 '+s_lgd_1)
    #         # axE.plot(dataN[lll_load][0,:st0_stop:2],eigs1[:st0_stop:2,0],color[2]+'_',label='eig0, no treatment')
    #         # axE.plot(dataN[lll_load][0,:st0_stop:2],eigs1[:st0_stop:2,1],color[2]+'|',label='eig1, no treatment')
    #         axE.legend(fancybox=True, framealpha=0.5);
    #         axE.set_ylabel("eigs per sim")
    #         axE.set_xlabel("time (days)")
    #         #axE.set_ylim(-.02,.001)
    #         # axE.set_yscale('log')
    #         axE.set_title(s_title)
    
    #for i, c in enumerate(cases):
    # data_drty = base_dirGD+"\\data\\"+base_model_name+"\\"+model_suffix+"\\"+case+"\\"+date_data_dir+ deathVal_dir + sub_drty;
    #G:\\My Drive
    #C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)
    
    # use_srvQ = True;
    # if use_srvQ:
    #     sig_vec = [[0,0],[3.6,0.05],[0.05, 3.6]];
    #     s_srv_vec = np.arange(len(sig_vec)).tolist();
    #     case = 'mu_test3\\26_Oct';
    #     for vv, el in enumerate(sig_vec):
    #         s_srv_vec[vv] = '\\' + str(el[0]) + '\\' + str(el[1]) + '\\';
    # else:
    #     sig_vec=[[3.6,0.05]];
    #     s_srv_vec = [''];
