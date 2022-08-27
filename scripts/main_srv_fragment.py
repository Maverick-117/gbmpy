# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 12:06:25 2022

@author: jhvo9
"""
import numpy as np
import funciones
from os import makedirs
from os.path import exists
#### SETUP

saveQ = False;
use_muQ = True;
subSelectQ = True; # true: downselect on a particular kind of feedback regime; false: use it all
compDosesQ = True; #true: yes you're comparing doses; false: no, you're comparing reprogramming versus non-reprogramming
kimICQ = True; # true: use kim's IC; false: use Yu's IC
c_dep_sQ = False; # true: make c in radiotherapy dependent on s; false: don't do that
kimReprogQ = False; 
kimDeathValQ = True; # true: sets death value to what Kim used; false: sets death value equal to 
    
model0Q = True; # true: using vo's no-survivin model and equations; false: not using it 
model1Q = False; # true: using vo's survivin-protects model and equations; false: not using it
model2Q = False; # true: using vo's survivin-protects-dedifferentiates model and equations; false: not using it

### VARS RELATED TO EXPERIMENTAL QUESTS ###
goldenLaptopQ = False;
schedChangeQ = False;
deathFdbkQ = False; # false: set death feedback gain to 0; true: don't use the false option
#useMuQ = True;

# Radiotherapy model
a,b =  np.array([0.17, 0.02]); 
a1 = 0.01
b1 = 1.77e-7
a2 = 0.125
b2 = 0.028
F = 0.016
#sig =  1e10;
DT =  3.9; # doubling time in days in order labeled by varible "cell_lines"
treat_start = 100; # treatment start time relative to start of simulation
acq_days_after_RT = 110; # simulation length after last fraction

# DE parameters
r1 = np.log(2)/DT; # growth rate of CSC
r2 = np.log(2)/DT; # growth rate of DCC
p = .505; # self renewal probability 
l_w = 10**(-7); # weak feedback on prob
l_s = 10**3; # strong feedback on prob
h1 = 10**5; # feedback on css div 
h2 = 10**5; # feedback on dcc div
pwr = 3;#Inf;
ROI_radius = 1; # Radius of the simulation region of intrest
rho = 10**9; # density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*np.pi*(ROI_radius ** 3)*rho;
post_therapy_end = 150000;
#mu_start = 0.0143;
ss = [1];#np.arange(12).tolist();#np.arange(12).tolist();#[6];#[4,8,9,5,10,11];
#[0,1,3,4,8,9,5,10,11];#[4,8,9,5,10,11];[0,1,2,3,6,7]
z = 1; n = 1;
l_vec =  [0,  0, l_w, l_s, l_w , l_s , 0 ,0 ,l_w,l_w,l_s,l_s];
h1_vec = [0, h1, 0  , 0  , h1  , h1  , h1,0 ,h1 ,0  , h1, 0 ];
h2_vec = [0, h2, 0  , 0  , h2  , h2  , 0 ,h2,0  ,h2 , 0 , h2];

# man-behind-the-curtain setup.
switch_vec = [subSelectQ, use_muQ, compDosesQ, deathFdbkQ, c_dep_sQ, kimReprogQ, kimDeathValQ, kimICQ, model0Q, model1Q, model2Q];
misc_pars = [DT, post_therapy_end, pwr, h1, h2, l_vec, ss, a, b];
par_setup_vec, string_setup_vec = funciones.parameter_setup(switch_vec, misc_pars);
total_start_frac, d, Doses, Frac, C, cont_p_a, rho, mu_bar, hd, cont_c, rng, mu_start, xi1, xi2, filler = par_setup_vec;
hd_str_mod, v_suffix, v_suffix_save_A, v_suffix_save, fg_subtitle, reprog_suffix, hd_str, hd_suffix, group_name, c_fdbk_mod = string_setup_vec;

total_start_cell_num = total_cell_num*total_start_frac;
cont_p_b = 10*cont_p_a; compt_mult = 100; # control the strength of inhibitory signal
srvn_zeta = [3.6, 0.05]; # defunct
srvn_csc = srvn_zeta[0]; srvn_dcc = srvn_zeta[1]; # defunct

if use_muQ:
    sig_vec = [0.1];#sorted(list(map(lambda p:  10 ** (p), np.arange(-1,3).tolist()))+[0]);
    rho_vec = [0.2];#sorted(list(map(lambda p: 2 ** (p), np.arange(0,-3,-1).tolist())) + [0,2]);
# rename rho to rho
    xi1_vec = [0.01];#sorted(list(map(lambda p:  10 ** (p), np.arange(-1,3).tolist()))+[0]);
    xi2_vec = [1];#sorted(list(map(lambda p:  10 ** (p), np.arange(-1,3).tolist()))+[0]);#np.arange(0.01,0.1,0.01);#[.1];#sorted(list(map(lambda p:  10 ** (p), np.arange(0,3).tolist()))+[0])
    xi3_vec = [1e9];
else:
    sig_vec = [0];

if schedChangeQ:
    Nsched_vec = [1];
else:
    Nsched_vec = [0];
# Nsched_vec = [0,1];
# stray strings to setup
fdbk_type_vec = ['No', 'Div Only', 'Weak', 'Strong', 'Weak', 'Strong'];
title_vec = [" nothing", " m_u, m_v", " p (weak)", " p (strong)", " m_u, m_v" + " and " + " p (weak)", " m_u, m_v" + " and " + " p (strong)"] + [hd_suffix];
color = ['k','r','b','m'];
cell_lines = ["U373MG"];

#### VARIATION ZONE FOR PARAMETERS
d *= 1;
time_pts1 = 500;
time_pts2 = 500;

day_month = "25_Feb"; #05_July last used
base_model_name = 'k2_model' # note: all the rho-sig sims should really fall into a model folder...
model_suffix = "conventional_BED\\"#"comparing_conventionals\\";
case = "reprogCheck";#"baseline_rho" or "mu_besi" 
if goldenLaptopQ:
    base_dirGD = "/DFS-L/DATA/lowengrub/vojh1";#"C:\\Users\\jhvo9\\Google Drive (vojh1@uci.edu)"
else:
    base_dirGD = "G:\\My Drive"
drty = base_dirGD + "\\a PhD Projects\\GBM Modeling\\python scripts\\data\\"+base_model_name+"\\"+model_suffix+"\\"+case+"\\"+day_month; # _div_rate_diff

if kimDeathValQ:
    deathVal_dir = '\\death_val_of_kim';
else:
    deathVal_dir = '\\death_val_of_not_kim';
if kimReprogQ:
    sub_drty = "\\kim_reprog";
else:
    sub_drty = "\\corrected_reprog";
total_drty = drty + deathVal_dir + sub_drty;
if not exists(total_drty):
    makedirs(total_drty);
if compDosesQ:
    comp_str = "Gy";
    comp_list = Doses;
    Frac_list = Frac;
else:
    comp_str = "reprog";
    comp_list = C;
    Frac_list = [Frac, Frac];
###
if deathFdbkQ:
    deathFdbk_str="_w_death_fdbk";
else:
    deathFdbk_str="_w_no_death_fdbk";
print(deathVal_dir + '\\' + sub_drty+comp_str+deathFdbk_str);
if compDosesQ:
    variegator = Doses; #outer loop will put the unchanging vector first
    fixor = C; # inner loop will use what's changing
    # technically this is unnecessary for data generation
    # but for plotting this helps
    # this could be easily generalized for any C and Doses used
else:
    variegator = C;
    fixor = Doses;
case_length = len(variegator);

for lll in rng:
    print("index val:", lll,"\n")
    u_sc = list(range(case_length)); u_dc = list(range(case_length)); u_srv = list(range(case_length)); t_vec = list(range(case_length));
    un_sc = list(range(case_length)); un_dc = list(range(case_length)); un_srv = list(range(case_length)); tn_vec = list(range(case_length));
    l = l_vec[lll]; h1 = h1_vec[lll]; h2 = h2_vec[lll];
    ## Tumor Growth ODE and radiotherapy simulation
    for useMuQ in [True, False]:
        s_useMuQ = 'w_reprog' if useMuQ else 'w_no_reprog';
        for ss, els in enumerate(xi3_vec):
            xi3 = els;
            s_xi3 = '\\'+str(els)+'\\';
            for tt, elt in enumerate(xi1_vec):
                xi1 = elt;
                s_xi1 = '\\'+str(elt)+'\\';
                for uu, elu in enumerate(xi2_vec):
                    xi2 = elu;
                    s_xi2 = '\\'+str(elu)+'\\';
                    for ww, elb in enumerate(rho_vec):
                        rho = elb;
                        s_rho = '\\' + str(elb) + '\\';
                        for vv, el in enumerate(sig_vec):
                            sig = el;
                            s_sig = '\\' + str(el) + '\\';
                            for nn, N_sched in enumerate(Nsched_vec):
                                total_drty_rho_sig = total_drty+s_xi1+s_xi2+s_xi3+ s_rho + s_sig+"\\Schedule"+str(N_sched);
                                if not exists(total_drty_rho_sig):
                                    makedirs(total_drty_rho_sig);
                                for gg in range(len(cell_lines)):
                                    print('Cell Line: '+cell_lines[gg]+'\n')
                                    print('a = '+str(a)+' b = '+str(b)+'\n');
                                    for g in range(len(fixor)):
                                        if compDosesQ:
                                            fix_str = str(C[g])+"_reprog";
                                        else:
                                            fix_str = str(Doses[g])+"_Gy";
                                            
                                        print("g = ",g)
                                        if compDosesQ:
                                            c = C[g];
                                        else:
                                            D = Doses[g];
                                            frac_num = Frac[g];
                                        for k in range(case_length):
                                            print("k = ",k)
                                            ## looping over reprogramming values
                                             # Defining Variables
                                            if compDosesQ:
                                                D = Doses[k];
                                                frac_num = int(Frac[k]);
                                            else:
                                                c = C[k];
                                                
                                            sc_start = total_start_frac*F;
                                            tc_start = total_start_frac-sc_start;
                                            total_cell_num = 4/3*np.pi*10**9
                                            total_start_num = total_start_frac*total_cell_num
                                            
                                            # Defining treatment days and ODE simulation start time after each fraction
                                            weeks = int(frac_num//5);
                                            total_days = frac_num + 2*weeks;
                                            acq_end = treat_start + post_therapy_end;#treat_start + total_days + acq_days_after_RT - 1;
                                            treat_days = np.tile([1,1,1,1,1,0,0], weeks+1)[:total_days].nonzero()[0]+treat_start;
                                            # if N_sched >= 1:
                                            #     if frac_num > 1:
                                            #         treat_days[1] = treat_days[-1]+1;
                                            #         treat_days.sort();
                                            # treat_days = np.tile([3,0,0,1,1,0,0], weeks+1)[:total_days].nonzero()[0]+treat_start;
                                            # for i in range(len(treat_days2)):
                                            #     if i % 3 == 0:
                                            #         sim_resume_days2.append();
                                            #     else:
                                            #         sim_resume_days2[i]+=10/(60*24);
                                            sim_resume_days = treat_days+10/(60*24); # ODE simulation resumes 10 minutes after RT
                                            treat_days = np.array(treat_days.tolist() + [acq_end]);  
                                            LQ_param = [a1, b1, a2, b2, c, D];
                                            surv_vec = [cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc, cont_c, useMuQ]; #assuming these control parameters are constant                                        
                                            sim_values = [model0Q, model1Q, model2Q, kimReprogQ, total_cell_num, treat_days, mu_start, LQ_param, total_start_frac, sc_start, sim_resume_days, surv_vec, time_pts1, time_pts2];
                                            para_values = (r1, r2, d, p, h1, h2, hd, z, l, n, sig, mu_bar, rho, xi1, xi2, D, xi3);
                                            
                                            U, T, U_none, T_none = funciones.dynamics(para_values,sim_values);
                                            u_sc[k] = U[0,:];
                                            u_dc[k] = U[1,:];
                                            u_srv[k] = U[2,:];
                                            t_vec[k] = T;
                                            un_sc[k] = U_none[0,:];
                                            un_dc[k] = U_none[1,:];
                                            un_srv[k] = U_none[2,:];
                                            tn_vec[k] = T_none;

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

load_temp = [1];#np.arange(12).tolist();#[1,3,5,6,7,10,11];#np.arange(12).tolist();


#[1,6,7] had...the same stability behaviour/eigenvalue plots
saveQ = False;
cases = ["2 Gy", "2.67 Gy" , "3.4 Gy", "5 Gy"]; 
comp_str = "Gy";
Doses = [10.0,10.0];
Fracs = [calc_BED_Frac(.17,.02,Doses[0]),calc_BED_Frac(.17,.02,Doses[1])]#list(map(lambda d:calc_BED_Frac(a,b,d,40*(1+(3.4)/8.5)),Doses));#;[25,20]
EOTs  = [calc_EOT(Fracs[0],100), calc_EOT(Fracs[1],100)];
EOT_diff = EOTs[0] - EOTs[1];
EOT_str = ",EOT-shifted"
match_EOT =True;
if not match_EOT:
    EOT_diff = 0.;
    EOT_str = "";
C = [5.196*10**(-3), 5.196*10**(-3)];
comp_dir = str(C[0])+"_reprog";
fix_str= str(Doses[-1])+"Gy"
comp_list = C;
unit = ' Gy';
T_stop = post_therapy_end;

log_yscaleQ = True;
if log_yscaleQ:
    log_str = " (log)";
    log_dir = "\\log_yscale\\";
else:
    log_str = "";
    log_dir = "\\linear_yscale\\"
#T_stop = 300;
T_stop = T_stop - 100;
total_cell_num = 4/3*np.pi*(1 ** 3)*1e9;
s_sched0 = str((Doses[0],Fracs[0]));
s_sched1 = str((Doses[1],Fracs[1]));
s_reprg0 = f"{float(C[0]):.3}";
s_reprg1 = f"{float(C[1]):.3}";
s_reprg_title = ";c:";

s_sched_title = ";sched:"
s_title_tweakable = s_sched_title+s_sched0;
s_lgd_0 = s_sched0;#str(csr[0])+unit+EOT_str;
s_lgd_1 = s_sched1+";no reprog;"+EOT_str;#str(csr[1])+unit;
COL2 = MplColorHelper('viridis', 0, len(Nsched_vec)-1)
COL1 = MplColorHelper('cool', 0, len(Nsched_vec)-1)
showNonReprogQ = False;

s_xi1_vec = ['\\' + str(el) for el in xi1_vec];
s_xi2_vec = ['\\' + str(elb) for elb in xi2_vec];
s_xi3_vec = ['\\' + str(els) for els in xi3_vec];
s_sig_vec = ['\\' + str(el) for el in sig_vec];
s_rho_vec = ['\\' + str(elb) + '\\' for elb in rho_vec];
load_temp = [1];
for lll_load in load_temp:
    h1 = h1_vec[lll_load];
    h2 = h2_vec[lll_load];
    l = l_vec[lll_load];
    figCellNo, axCellNo = plt.subplots(ncols=2, nrows=2, figsize=(20, 20),
                        constrained_layout=True)
    figRest, axRest = plt.subplots(ncols=2, nrows=2, figsize=(20, 20),
                            constrained_layout=True)
    #figCellNo, axCellNo[0,0] = plt.subplots(figsize=(8,8));
    #figCellNo, axCellNo[0,1] = plt.subplots(figsize=(8,8));
    s_title = "(l:"+f"{float(l_vec[lll_load]):.3}"+";h1:"+f"{float(h1_vec[lll_load]):.3}"+";h2:"+f"{float(h2_vec[lll_load]):.3}"+";hd:"+f"{float(hd):.3}"+s_title_tweakable+")";
    # s_lgd_0 = "Schedule "+str(eln)+Nsched_vec_str[nn]+';'+str((Doses[nn],Fracs[nn]));
    # s_lgd_1 = "w/ constant reprog";#"Schedule "+str(eln)+Nsched_vec_str[nn]+'alt';

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
                            st2_stop = np.flatnonzero(T < T_stop).max();
                            axCellNo[0,0].plot(T[:st2_stop],U[0,:st2_stop] * total_cell_num,'--',color=COL2.get_rgb(nn),label='CSC '+s_lgd_0+'_'+str(lll_load))
                            axCellNo[0,1].plot(T[:st2_stop],U[1,:st2_stop] * total_cell_num,'-',color=COL2.get_rgb(nn),label='DCC '+s_lgd_0+s_sig+'_'+str(lll_load))
                            if showNonReprogQ:
                                st1_stop = np.flatnonzero(T_none < T_stop).max();
                                axCellNo[0,0].plot(T_none[:st2_stop]+EOT_diff,U_none[0,:st2_stop] * total_cell_num,'--',color=COL1.get_rgb(nn),label='CSC '+s_lgd_1+s_sig+'_'+str(lll_load))
                                axCellNo[0,1].plot(T_none[:st1_stop]+EOT_diff,U_none[1,:st1_stop] * total_cell_num,'-',color=COL1.get_rgb(nn),label='DCC '+s_lgd_1+s_sig+'_'+str(lll_load))
                    
    st0_stop = np.flatnonzero(T_none < T_stop).max();
    axCellNo[0,0].plot(T_none[:st0_stop],U_none[0,:st0_stop] * total_cell_num,'g--',label='CSC, no treatment_'+str(lll_load))
    axCellNo[0,1].plot(T_none[:st0_stop],U_none[1,:st0_stop] * total_cell_num,'g-',label='DCC, no treatment_'+str(lll_load))
    axCellNo[0,0].legend(fancybox=True, framealpha=0.5);
    axCellNo[0,0].set_ylabel("CSC Cell Number"+log_str)
    axCellNo[0,0].set_xlabel("time (days)")
    axCellNo[0,0].set_title(s_title)
    axCellNo[0,1].legend(fancybox=True, framealpha=0.5);
    axCellNo[0,1].set_ylabel("DCC Cell Number"+log_str)
    axCellNo[0,1].set_xlabel("time (days)")
    axCellNo[0,1].set_title(s_title)
    if log_yscaleQ:
        axCellNo[0,0].set_yscale('log')
        axCellNo[0,1].set_yscale('log')
    if saveQ:
        figCellNo.savefig(plot_drty+"pop"+plot_save_suffix+".png",dpi=300)
            #replace("pop"+plot_save_suffix+".png",plot_drty+"pop"+plot_save_suffix+".png")

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
                        
                        st2_stop = np.flatnonzero(T < T_stop).max();
                        
                        # if log_yscaleQ:
                        #     axRest[0,0].set_yscale('log')
                        axRest[0,0].plot(T[:st2_stop],U[2,:st2_stop],'--',color=COL2.get_rgb(nn),label='mu frac'+s_lgd_0+s_sig+'_'+str(lll_load))
                        if showNonReprogQ:
                            st1_stop = np.flatnonzero(T_none < T_stop).max();
                            axRest[0,0].plot(T_none[:st1_stop]+EOT_diff,data1[nn][tt][uu][ss][ww][vv][lll_load][3,:st1_stop],'--',color=COL1.get_rgb(nn),label='mu frac'+s_lgd_1+s_sig+'_'+str(lll_load))

    st0_stop = np.flatnonzero(T_none < T_stop).max();
    axRest[0,0].plot(T_none[:st0_stop],U_none[2,:st0_stop],'g--',label='mu frac, no treatment_'+str(lll_load))
    axRest[0,0].legend(fancybox=True, framealpha=0.5);
            
    axRest[0,0].set_ylabel("mu frac"+log_str)
    axRest[0,0].set_xlabel("time (days)")
    axRest[0,0].set_title(s_title)
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
                        
                        st2_stop = np.flatnonzero(T < T_stop).max();
                        axRest[0,1].plot(T[:st2_stop],beta * U[2,:st2_stop] * U[1,:st2_stop] * total_cell_num,'--',color=COL2.get_rgb(nn),label='mu rate'+s_lgd_0+s_sig+'_'+str(lll_load))
                        axRest[0,1].set_ylabel("mu rate"+log_str)
                        axRest[0,1].set_xlabel("time (days)")
                        axRest[0,1].set_title(s_title)
                        if showNonReprogQ:
                            st1_stop = np.flatnonzero(T_none < T_stop).max();
                            axRest[0,1].plot(T_none[:st1_stop]+EOT_diff,beta * data1[nn][tt][uu][ss][ww][vv][lll_load][3,:st1_stop] * U_none[1,:st1_stop] * total_cell_num,'--',color=COL1.get_rgb(nn),label='mu rate'+s_lgd_1+s_sig+'_'+str(lll_load))
    
    st0_stop = np.flatnonzero(T_none < T_stop).max();
    axRest[0,1].plot(T_none[:st0_stop],beta * U_none[2,:st0_stop] * U_none[1,:st0_stop] * total_cell_num,'g--',label='mu rate, no treatment_'+str(lll_load))
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
                        
                        st2_stop = np.flatnonzero(T < T_stop).max();
                        axRest[1,0].plot(T[:st2_stop],xi3* U[0,:st2_stop]/(1+xi3* U[0,:st2_stop]),'--',color=COL2.get_rgb(nn),label='feedforward rate'+s_lgd_0+s_sig+'_'+str(lll_load))
                        axRest[1,0].set_ylabel("feedforward rate"+log_str)
                        axRest[1,0].set_xlabel("time (days)")
                        axRest[1,0].set_title(s_title)
                        if showNonReprogQ:
                            st1_stop = np.flatnonzero(T_none < T_stop).max();
                            axRest[1,0].plot(T_none[:st1_stop]+EOT_diff,xi3 * U_none[0,:st1_stop]/(1+ xi3*U_none[0,:st1_stop]),'--',color=COL1.get_rgb(nn),label='feedforward rate'+s_lgd_1+s_sig+'_'+str(lll_load))
    
    st0_stop = np.flatnonzero(T_none < T_stop).max();
    axRest[1,0].plot(T_none[:st0_stop],xi3 * U_none[0,:st0_stop]/(1+xi3*U_none[0,:st0_stop]),'g--',label='feedforward rate, no treatment_'+str(lll_load))
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
                        
                        st2_stop = np.flatnonzero(T < T_stop).max();
                        #
                        axRest[1,1].plot(T[:st2_stop],r2/(1+h2*U[1,:st2_stop]) * xi3*U[0,:st2_stop]/(1+xi3*U[0,:st2_stop]) - d-(1.1*r2-d)*hd*U[1,:st2_stop]/(1+hd*U[1,:st2_stop]),'--',color=COL2.get_rgb(nn),label='effective division rate'+s_lgd_0+s_sig+'_'+str(lll_load))
                        axRest[1,1].set_ylabel("effective division rate"+log_str)
                        axRest[1,1].set_xlabel("time (days)")
                        axRest[1,1].set_title(s_title)
                        if showNonReprogQ:
                            st1_stop = np.flatnonzero(T_none < T_stop).max();
                            axRest[1,1].plot(T_none[:st1_stop]+EOT_diff,r2/(1+h2*U_none[1,:st1_stop]) * xi3*U_none[0,:st1_stop]/(1+xi3*U_none[0,:st1_stop]) - d-(1.1*r2-d)*hd*U_none[1,:st1_stop]/(1+hd*U_none[1,:st1_stop]),'--',color=COL1.get_rgb(nn),label='effective division rate'+s_lgd_1+s_sig+'_'+str(lll_load))
    
    st0_stop = np.flatnonzero(T_none < T_stop).max();
    axRest[1,1].plot(T_none[:st0_stop],r2/(1+h2*U_none[1,:st0_stop]) * xi3*U_none[0,:st0_stop]/(1+xi3*U_none[0,:st0_stop]) - d-(1.1*r2-d)*hd*U_none[1,:st0_stop]/(1+hd*U_none[1,:st0_stop]),'g--',label='effective division rate, no treatment_'+str(lll_load))
    axRest[1,1].legend(fancybox=True, framealpha=0.5);


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
                        
                        st2_stop = np.flatnonzero(T < T_stop).max();
                        
                        if log_yscaleQ:
                            axCellNo[1,0].set_yscale('log')
                        axCellNo[1,0].plot(T[:st2_stop],(U[0,:st2_stop]+U[1,:st2_stop]) * total_cell_num,'--',color=COL2.get_rgb(nn),label='total '+s_lgd_0+s_sig+'_'+str(lll_load))
                        axCellNo[1,0].set_ylabel("Total Cell Number"+log_str)
                        axCellNo[1,0].set_xlabel("time (days)")
                        axCellNo[1,0].set_title(s_title)
                        if showNonReprogQ:
                            st1_stop = np.flatnonzero(T_none < T_stop).max();
                            axCellNo[1,0].plot(T_none[:st1_stop]+EOT_diff,(U_none[0,:st1_stop]+U_none[1,:st1_stop]) * total_cell_num,'--',color=COL1.get_rgb(nn),label='total '+s_lgd_1+s_sig+'_'+str(lll_load))

    st0_stop = np.flatnonzero(T_none < T_stop).max();
    axCellNo[1,0].plot(T_none[:st0_stop],(U_none[0,:st0_stop]+U_none[1,:st0_stop]) * total_cell_num,'g--',label='total, no treatment_'+str(lll_load))
    axCellNo[1,0].legend(fancybox=True, framealpha=0.5);
    # if saveQ:
    #     figCellNo.savefig("total_pop"+plot_save_suffix+".png",dpi=300)
    #     # replace(base_dirGD+"total_pop"+plot_save_suffix+".png",plot_drty+"total_pop"+plot_save_suffix+".png")
        
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
                        csc_p += [U[0]/(U[0] + U[1])];
                        csc_p += [U_none[0]/(U_none[0] + U_none[1])];
                        
                        st2_stop = np.flatnonzero(T < T_stop).max();
                        
                        
                        # if log_yscaleQ:
                        #     axCellNo[1,1].set_yscale('log')
                        axCellNo[1,1].plot(T[:st2_stop],csc_p[0][:st2_stop],'-',color=COL2.get_rgb(nn),label='CSC frac'+s_lgd_0+s_sig+'_'+str(lll_load))
        
                        if showNonReprogQ:
                            st1_stop = np.flatnonzero(T_none < T_stop).max();
                            axCellNo[1,1].plot(T_none[:st1_stop]+EOT_diff,csc_p[1][:st1_stop],'-',color=COL1.get_rgb(nn),label='CSC frac'+s_lgd_1+s_sig+'_'+str(lll_load))
        # replace(base_dirGD+"csc"+plot_save_suffix+".png",plot_drty+"csc"+plot_save_suffix+".png")
    csc_p += [U_none[0]/(U_none[0] + U_none[1])];
    st0_stop = np.flatnonzero(T_none < T_stop).max();
    axCellNo[1,1].plot(T_none[:st0_stop],csc_p[2][:st0_stop],'g-',label='CSC frac, no treatment')
    axCellNo[1,1].legend(fancybox=True, framealpha=0.5);
    axCellNo[1,1].set_ylabel("CSC frac")
    axCellNo[1,1].set_xlabel("time (days)")
    axCellNo[1,1].set_title(s_title)
    # if saveQ:
    #     figC.savefig("csc"+plot_save_suffix+".png",dpi=300)
