# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 20:11:28 2023

@author: jhvo9
"""

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

def dU_dt_kim(U,t, r1, r2, d, p, h1, h2, hd, z, l, n):
    return np.array([(2*p/(1+l*U[1]**n)-1)*r1/(1+h1*U[1]**n)*U[0], # stem cell
                     2*(1-p/(1+l*U[1]**n))*r1/(1+h1*U[1]**n)*U[0]+r2/(1+h2*U[1]**z)*U[1]-d*U[1],
                     0
                     ])
def calc_BED_Frac(a,b,Dose,baseline_BED=60*(1+2/8.5)):
    # a = .17, b =.02
    # the baseline BED is set to conventional dosage, by default
    if Dose > 0:
        return np.floor(baseline_BED/(Dose*(1+Dose/(a/b))))
    else:
        return 0;
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
C = np.array([0, 5.196*10**(-3)]);

r1 = np.log(2)/DT; # growth rate of CSC
r2 = np.log(2)/DT; # growth rate of DCC
d = r1/5;
p = .505; # self renewal probability 
l_w = 10**(-7); # weak feedback on prob
l_s = 10**3; # strong feedback on prob
h1 = 0; # feedback on css div 
h2 = 0; # feedback on dcc div
l = 1e3;
hd = 0; z = 1; n = 1;



total_start_frac = 0.0005/64;
ROI_radius = 1; # Radius of the simulation region of intrest
rho = 10 ** 9; # density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*np.pi*ROI_radius**3*rho;
total_start_cell_num = total_cell_num*total_start_frac;
Doses = [2.0];#np.arange(1,12,dtype=float);#np.arange(2,22,2,dtype=float);#sorted(np.arange(1,21,dtype=float).tolist() + [40/15,34/10]);
Frac = list(map(lambda d: float(np.floor(calc_BED_Frac(.17,.02,d))),Doses)); # 25*(1+(5)/8.5)
post_therapy_end = 200;
for c in C:        
    sc_start = total_start_frac*F;
    tc_start = total_start_frac-sc_start;
    total_cell_num = 4/3*np.pi*10**9
    total_start_num = total_start_frac*total_cell_num
    
    # Defining treatment days and ODE simulation start time after each fraction
    frac_num = Frac[0];
    D = Doses[0];
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
    para_values = (r1, r2, d, p, h1, h2, hd, z, l, n);
    U0 = [sc_start, tc_start,0]; 
    T = np.linspace(0, post_therapy_end, 1000);
    U = integrate.odeint(dU_dt_kim, U0, T, rtol=1e-10, atol=1e-10, args=para_values).T
    
    plt.plot(T,U[0]*total_cell_num,'g--');
    plt.plot(T,U[1]*total_cell_num,'g-');
    plt.yscale('log')
