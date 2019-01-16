# -*- coding: utf-8 -*-

'''
# MATLAB Code written by: Pierre Desir
# Adapted in Python by: Yifan Wang
# Date: 01-16-2019
 
#This MatLab Code solves the PFR model for the acid-catalyzed dehydration of
#fructose to HMF using HCl as the catalyst and generates the kinetic data
#for conversion, yield, and selectivity of the species as a function of
#temperature, pH, and time. All the kinetic and thermodynamic parameters 
#were taken from the kinetic model by T. Dallas Swift et al. ACS Catal 
#2014, 4,259-267.
'''
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

#%%------------------------------------------------------------------------
# INPUT PARAMETERS: n_T, pH, Tmin_degC, Tmax_degC, tmin, tmax, Fru0
n_T = 11 #Number of temperature points
pH = 0.7 #Rxn pH
Tmin_degC = 100  #Min rxn temperature [°C]
Tmax_degC = 200  #Max rxn temperature [°C]
t0 = 0 #Initial time point [min]
tf = 1e3 #Final time point[min]
Fru0 = 1 #Normalized initial fructose concentration (always equal to 1)
#-------------------------------------------------------------------------
#VARIABLE REACTION PARAMETERS
C_Hplus = 10**(-pH) #H+ concentraction [mol/L]
T_degC = np.linspace(Tmin_degC,Tmax_degC,n_T) #Rxn temperature [°C]
T_K = T_degC + 273*np.ones(len(T_degC)) #Rxn temperature [K]
#-------------------------------------------------------------------------
#%%
def PFR(C, t, T_K):
    '''
    #The "PFR" function describes the species conservative equation as a
    #set of ODEs at T_K
    '''    
    #CONSTANTS
    R = 8.314 #Ideal gas law constant [J/mol-K]
    Tm = 381 #Mean temperature of all Rxn [K]
    
    #OTHER MODEL PARAMETERS
    C_H2O = 47.522423477254065 + 0.06931572301966918*T_K 
    - 0.00014440077466393135*T_K**2 #Water 
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
    K_bp = K_bp303*np.exp(-(delH_bp/R)*(1/T_K-1/303))
    K_bf = K_bf303*np.exp(-(delH_bf/R)*(1/T_K-1/303))
    K_ap = K_ap303*np.exp(-(delH_ap/R)*(1/T_K-1/303))
    K_af = K_af303*np.exp(-(delH_af/R)*(1/T_K-1/303))
    
    #Furanose fraction at equilibirum as a function of temperature
    phi_f = (K_af + K_bf/(1 + K_af + K_bf + K_ap + K_bp))
    
    # ACTIVATIONS ENERGIES FOR RXN1,RXN2,...,RXN5
    Ea = np.array([127, 133, 97, 64, 129]) * (10**3) #[J/mol]
    
    #NTURAL LOG OF RXN RATE CONSTANTS AT 381 K FOR RXN1,RXN2,...,RXN5
    lnk381 = np.array([1.44, -4.22, -3.25, -5.14, -4.92])
     
    #RXN RATE CONSTANTS FOR RXN1,RXN2,...,RXN5 AS A FUNCTION OF 
    #TEMPERATURE
    k = np.zeros(5)
    for i in range(len(k)):
        k[i] = np.exp(lnk381[i]-(Ea[i]/R)*(1/T_K-1/Tm)) #[min^-1.M^-1]
    
    #RXN RATES FOR THE RXN NETWORK OF FRUCTOSE DEHYDRATION
    #Note: C[0] = Normalized Fructose concentration; C[1] = Normalized
    #HMF concentration; C[2] = Normalized LA concentration; 
    #C[3] = Normalized FA concentration;
    Rxn = np.zeros(5) #[mol/L-min]
    Rxn[0] = k[0]*phi_f[0]*C[0]*C_Hplus/C_H2O #[mol/L-min]
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
  

#def model(T):
    
    # Solving for the PFR model at certain temperature