# -*- coding: utf-8 -*-

'''
# MATLAB Code written by: Pierre Desir
# Adapted in Python by: Yifan Wang
# Date: 01-16-2019

# Plots of reactor performance at 11 temperatures are shown in this code

 
#This MatLab Code solves the PFR model for the acid-catalyzed dehydration of
#fructose to HMF using HCl as the catalyst and generates the kinetic data
#for conversion, yield, and selectivity of the species as a function of
#temperature, pH, and time. All the kinetic and thermodynamic parameters 
#were taken from the kinetic model by T. Dallas Swift et al. ACS Catal 
#2014, 4,259-267.
'''
import time 
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.integrate import odeint, ode

t = time.time()
#%%-------------------------------------------------------------------------
#Reactor functions
def Reactor(**conditions):
    
    '''
    Reactor model under certain temperature, pH and residence time tf
    '''
    T_degC = conditions['T_degC'] # reaction temperature in C
    pH = conditions['pH'] # reaction pH
    tf = conditions['tf'] #Final time point [min]
    t0 = 0 #Initial time point [min]
    Fru0 = 1 #Normalized initial fructose concentration (always equal to 1)

    def PFR(t,C, T_degC, pH):
        '''
        #The "PFR" function describes the species conservative equation as a
        #set of ODEs at T
        '''    
        #CONSTANTS
        R = 8.314 #Ideal gas law constant [J/mol-K]
        Tm = 381 #Mean temperature of all Rxn [K]
        

        #OTHER MODEL PARAMETERS
        T = T_degC + 273 # reaction temperarure in K
        C_Hplus = 10**(-pH) #H+ concentraction [mol/L]
        C_H2O = 47.522423477254065 + 0.06931572301966918*T - 0.00014440077466393135* (T**2) #Water 
        
        #concentration as a function of temperature from 25 °C to 300 °C
        #[mol/L]
           
        
        #TAUTOMER PARAMETERS
        #Enthalpy of Rxn b/t open-chain fructose and tautomers
        #Note: a = alpha; b = beta; p = pyrannose; f = furannose
        delH_bp = -30.2e3 #[J/mol]
        delH_bf = -19e3 #[J/mol]
        delH_ap = -5.5e3 #[J/mol]
        delH_af = -14.2e3 #[J/mol]
        
        #Equilibrium constants b/t open-chain fructose and tautomers at 303 K
        K_bp303 = 59.2
        K_bf303 = 26.4
        K_ap303 = 0.6
        K_af303 = 6.4
        
        #Equilibirium constants b/t open-chain fructose and tautomers as a 
        #function of temperature
        K_bp = K_bp303*np.exp(-(delH_bp/R)*(1/T-1/303))
        K_bf = K_bf303*np.exp(-(delH_bf/R)*(1/T-1/303))
        K_ap = K_ap303*np.exp(-(delH_ap/R)*(1/T-1/303))
        K_af = K_af303*np.exp(-(delH_af/R)*(1/T-1/303))
        
        #Furanose fraction at equilibirum as a function of temperature
        phi_f = (K_af + K_bf)/(1 + K_af + K_bf + K_ap + K_bp)
        
        # ACTIVATIONS ENERGIES FOR RXN1,RXN2,...,RXN5
        Ea = np.array([127, 133, 97, 64, 129]) * (10**3) #[J/mol]
        
        #NTURAL LOG OF RXN RATE CONSTANTS AT 381 K FOR RXN1,RXN2,...,RXN5
        lnk381 = np.array([1.44, -4.22, -3.25, -5.14, -4.92])
         
        #RXN RATE CONSTANTS FOR RXN1,RXN2,...,RXN5 AS A FUNCTION OF 
        #TEMPERATURE
        k = np.zeros(5)
        for i in range(len(k)):
            k[i] = np.exp(lnk381[i]-(Ea[i]/R)*(1/T-1/Tm)) #[min^-1.M^-1]
        
        #RXN RATES FOR THE RXN NETWORK OF FRUCTOSE DEHYDRATION
        #Note: C[0] = Normalized Fructose concentration; C[1] = Normalized
        #HMF concentration; C[2] = Normalized LA concentration; 
        #C[3] = Normalized FA concentration;
        Rxn = np.zeros(5) #[mol/L-min]
        Rxn[0] = k[0]*phi_f*C[0]*C_Hplus/C_H2O #[mol/L-min]
        Rxn[1] = k[1]*C[0]*C_Hplus #[min^-1]
        Rxn[2] = k[2]*C[1]*C_Hplus #[min^-1]
        Rxn[3] = k[3]*C[1]*C_Hplus #[min^-1]
        Rxn[4] = k[4]*C[0]*C_Hplus #[min^-1]
         
        #SPECIES CONSERVATIVE EQUATIONS
        #Notation: rhs = dC/dt
        rhs = np.zeros(4) #[mol/L-min]
        rhs[0] = (-Rxn[0]-Rxn[1]-Rxn[4]) #Fructose
        rhs[1] = (Rxn[0]-Rxn[2]-Rxn[3]) #HMF
        rhs[2] = Rxn[2] #LA
        rhs[3] = (Rxn[2]+Rxn[4]) #FA
        
        return rhs
      
    #%% SOLVING FOR THE PFR MODEL at certain temperature T_K
    C0 = np.array([Fru0, 0, 0, 0])
    
    # Construct the ode solver, ode45 with varying step size
    sol = []
    def solout(t, y):
        sol.append([t, *y])
    solver = ode(PFR).set_integrator('dopri5', rtol  = 1e-6, method='bdf')
    solver.set_solout(solout)
    #feed in argumenrs and initial conditions for odes
    solver.set_initial_value(C0, t0).set_f_params(*(T_degC, pH)) 
    solver.integrate(tf)
    sol = np.array(sol)
    
    Tau = sol[:,0] #time steps
    Conc = sol[:,1:] #concentration matrix
     
    #RESULTS
    Fru = np.around(Conc[:,0], decimals = 4) #Fructose normalized concentration 
    HMF = np.around(Conc[:,1], decimals = 4) #HMF normalized concentration 
    LA = np.around(Conc[:,2], decimals = 4) #Levulinic acid (LA) concentration 
    FA = np.around(Conc[:,3], decimals = 4) #Formic acid (FA) concentration 
    Conv = np.around(100*(1-Fru), decimals = 4) #Fructose conversion [%]
    HMF_Yield = 100*HMF #HMF yield [%]
    HMF_Select = 100*HMF_Yield/Conv #HMF selectivity [%]
    LA_Yield = 100*LA #LA yield [%]
    FA_Yield = 100*FA #FA yield [%]
     
    #OPTIMAL CONDITIONS FOR MAX HMF YIELD
    Max_HMF_Yield = max(HMF_Yield) #Maximum HMF Yield [%]
    index_range = np.where(HMF_Yield == Max_HMF_Yield)[0] # Index of matrix 
    #element where the HMF yield is at its max value
    index = index_range[max(np.where(Conv[index_range] == max(Conv[index_range]))[0])]
    #index of matrix element for optimal conditions for max HMF yield
    Tau_opt = Tau[index] #Optimal residence time to reach maximum HMF 
    #yield [min]
    Opt_Conv = np.around(Conv[index], decimals = 0) #Fructose conversion at max HMF yield [%]
    Opt_Select = HMF_Select[index] #HMF selectivity at max HMF yield [%]
    
    #REPORTING OPTIMAL CONDITIONS
    Opt_Cond = (T_degC, Tau_opt, Max_HMF_Yield, Opt_Conv, Opt_Select)
    #Temperature, optimal residence time, max HMF yield, conversion at max
    #HMF yield, HMF selectivity at max HMF yield


    return Tau, Conv, HMF_Select, HMF_Yield, LA_Yield, FA_Yield, Opt_Cond


#%%------------------------------------------------------------------------
# INPUT PARAMETERS: n_T, pH, Tmin_degC, Tmax_degC, tmin, tmax, Fru0
pH = 0.7 #Rxn pH
t0= 0  #Initial time point [min]
tf = 1e3 #Final time point[min]
Tmin_degC = 100 # %Min rxn temperature [°C]
Tmax_degC = 200 # %Max rxn temperature [°C]
n_T = 11; #Number of temperature points
#VARIABLE REACTION PARAMETERS
T_degC = np.linspace(Tmin_degC,Tmax_degC,n_T) #Rxn temperature [°C]
#-------------------------------------------------------------------------
# Run the reactor model and collect data at each temperature 
Tau = []
Conv = []
HMF_Select = []
HMF_Yield = []
LA_Yield = []
FA_Yield = []
Opt_Cond = []

for Ti in T_degC:
    Conditions = {'T_degC': Ti, 'pH': pH, 'tf' : tf}
    Tau_i, Conv_i, HMF_Select_i, HMF_Yield_i, LA_Yield_i, FA_Yield_i, Opt_Cond_i = Reactor(**Conditions)
    Tau.append(Tau_i)
    Conv.append(Conv_i)
    HMF_Select.append(HMF_Select_i)
    HMF_Yield.append(HMF_Yield_i)
    LA_Yield.append(LA_Yield_i)
    FA_Yield.append(FA_Yield_i)    
    Opt_Cond.append(Opt_Cond_i)
    
Opt_Cond = np.array(Opt_Cond)
#%% Plotting Section
# Set the parameters for plotting
font = {'size'   : 12}

matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.2
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['xtick.major.width'] = 1.2
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['ytick.major.width'] = 1.2

#FIG.1 : HMF selectivity vs. FRU conversion
fig = plt.figure(figsize=(8,6))
for i in range(n_T):
    plt.plot(Conv[i], HMF_Select[i], label = str(np.int(T_degC[i])) + r'$^oC$');
plt.xlabel('Fructose Conversion (%)')
plt.ylabel('HMF Selectity (%)')
plt.xlim(5,100)
plt.legend(loc = 'best', frameon=True)
plt.title('HMF Selectivity')
fig.savefig('HMF Selectivity'+'.png')

#FIG.2 : FRU Conversion vs. residence time
fig = plt.figure(figsize=(8,6))
for i in range(n_T):
    plt.semilogx(Tau[i], Conv[i], label = str(np.int(T_degC[i])) + r'$^oC$')
plt.xlabel('Residence Time (min)')
plt.ylabel('Fructose Conversion (%)')
plt.xlim(1e-4,tf)
plt.legend(loc = 'best', frameon=True)
plt.title('Fructose Conversion')
fig.savefig('Fructose Conversion'+'.png')

#FIG.3 : HMF yield vs. residence time
fig = plt.figure(figsize=(8,6))
for i in range(n_T):
    plt.semilogx(Tau[i], HMF_Yield[i], label = str(np.int(T_degC[i])) + r'$^oC$')
plt.xlabel('Residence Time (min)')
plt.ylabel('HMF Yield (%)')
plt.xlim(1e-4,tf)
plt.legend(loc = 'best', frameon=True)
plt.title('HMF Yield')
fig.savefig('HMF Yield'+'.png')

#FIG.4 : LA yield vs. residence time
fig = plt.figure(figsize=(8,6))
for i in range(n_T):
    plt.semilogx(Tau[i], LA_Yield[i], label = str(np.int(T_degC[i])) + r'$^oC$')
plt.xlabel('Residence Time (min)')
plt.ylabel('LA Yield (%)')
plt.xlim(1e-4,tf)
plt.legend(loc = 'best', frameon=True)
plt.title('LA Yield')
fig.savefig('LA Yield'+'.png')
    
#FIG.5 : FA yield vs. residence time
fig = plt.figure(figsize=(8,6))
for i in range(n_T):
    plt.semilogx(Tau[i], LA_Yield[i], label = str(np.int(T_degC[i])) + r'$^oC$')
plt.xlabel('Residence Time (min)')
plt.ylabel('FA Yield (%)')
plt.xlim(1e-4,tf)
plt.legend(loc = 'best', frameon=True)
plt.title('FA Yield')
fig.savefig('FA Yield'+'.png')

#FIG.6 : Max HMF yield vs. temperature
Max_HMF_Yield = Opt_Cond[:,2]
fig = plt.figure(figsize=(8,6))
plt.plot(T_degC, Max_HMF_Yield)
plt.xlabel('Temperature (°C)')
plt.ylabel('Max HMF Yield (%)')
plt.title('Max HMF Yield')
fig.savefig('Max HMF Yield'+'.png')

#FIG.7 : Optimal residence time vs. temperature
Tau_opt = Opt_Cond[:,1]
fig = plt.figure(figsize=(8,6))
plt.semilogy(T_degC, Tau_opt)
plt.xlabel('Temperature (°C)')
plt.ylabel('Optimal Residence Time (min)')
plt.title('Optimal Residence Time')
fig.savefig('Optimal Residence Time'+'.png')

#%%Obtain the CPU time
elapsed = time.time() - t
print('Elapsed time is %.5f seconds.'% elapsed)