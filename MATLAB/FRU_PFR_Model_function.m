function [Opt_Cond] = FRU_PFR_Model_function(T_degC, pH, tf)
% This code takes in 3 inputs
 
% T - Rxn temperature [°C]
% pH - Rxn pH 
% tf - final residence time (min)
% Example: FRU_PFR_Model_function(100, 0.7, 1000)


%Code written by: Pierre Desir
%Date: 06-13-16
 
%This MatLab Code solves the PFR model for the acid-catalyzed dehydration of
%fructose to HMF using HCl as the catalyst and generates the kinetic data
%for conversion, yield, and selectivity of the species as a function of
%temperature, pH, and time. All the kinetic and thermodynamic parameters 
%were taken from the kinetic model by T. Dallas Swift et al. ACS Catal 
%2014, 4,259-267.

 
 tic
 %-------------------------------------------------------------------------
 %INPUT PARAMETERS: n_T, pH, Tmin_degC, Tmax_degC, tmin, tmax, Fru0
 Tmin_degC = 100; %Min rxn temperature [°C]
 Tmax_degC = 200; %Max rxn temperature [°C]
 t0 = 0; %Initial time point [min]
 Fru0 = 1; %Normalized initial fructose concentration (always equal to 1)
 %-------------------------------------------------------------------------
 %VARIABLE REACTION PARAMETERS
 C_Hplus = 10^(-pH); %H+ concentraction [mol/L]
 T_K = T_degC + 273; %Rxn temperature [K]
 %-------------------------------------------------------------------------

    
     
 %SOLVING FOR THE PFR MODEL     
 options = odeset('RelTol',1e-6,'bdf','on' );
 [Tau,Conc]= ode15s(@PFR,[t0,tf],[Fru0 0 0 0 ],options);
 %Tau is the residence (space) time in the PFR [min]
 %Conc is the species concentration matrix [mol/L]

%RESULTS
Fru = round(Conc(:,1),4); %Fructose normalized concentration 
HMF = round(Conc(:,2),4); %HMF normalized concentration 
LA = round(Conc(:,3),4); %Levulinic acid (LA) concentration 
FA = round(Conc(:,4),4); %Formic acid (FA) concentration 
Conv = round(100*(1-Fru),4); %Fructose conversion [%]
HMF_Yield = 100*HMF; %HMF yield [%]
HMF_Select = 100*HMF_Yield./Conv; %HMF selectivity [%]
LA_Yield = 100*LA; %LA yield [%]
FA_Yield = 100*FA; %FA yield [%]

plot(Tau, HMF_Yield);

%OPTIMAL CONDITIONS FOR MAX HMF YIELD
Max_HMF_Yield = max(HMF_Yield); %Maximum HMF Yield [%]
index_range = (find(HMF_Yield == Max_HMF_Yield)); % Index of matrix 
%element where the HMF yield is at its max value
index = index_range(max(find(Conv(index_range) == max(Conv(index_range)))));
%index of matrix element for optimal conditions for max HMF yield
Tau_opt = Tau(index); %Optimal residence time to reach maximum HMF 
%yield [min]
Opt_Conv = round(Conv(index)); %Fructose conversion at max HMF yield [%]
Opt_Select = HMF_Select(index); %HMF selectivity at max HMF yield [%]

%REPORTING OPTIMAL CONDITIONS
Opt_Cond = [T_degC Tau_opt Max_HMF_Yield Opt_Conv Opt_Select];
%Temperature, optimal residence time, max HMF yield, conversion at max
%HMF yield, HMF selectivity at max HMF yield
    
 
%  
% toc
  
    function rhs = PFR(t,C)
        %The "PFR" function describes the species conservative equation as a
        %set of ODEs
        
        %CONSTANTS
        R = 8.314; %Ideal gas law constant [J/mol-K]
        Tm = 381; %Mean temperature of all Rxn [K]

        %OTHER MODEL PARAMETERS
        C_H2O = 47.522423477254065 + 0.06931572301966918*T_K...
           -0.00014440077466393135*T_K^2; %Water 
        %concentration as a function of temperature from 25 °C to 300 °C
        %[mol/L]
       
        
        %TAUTOMER PARAMETERS
        %Enthalpy of Rxn b/t open-chain fructose and tautomers
        %Note: a = alpha; b = beta; p = pyrannose; f = furannose
        delH_bp = -30.2e3; %[J/mol]
        delH_bf = -19e3; %[J/mol]
        delH_ap = -5.5e3; %J/mol]
        delH_af = -14.2e3; %[J/mol]
        
        %Equilibrium constants b/t open-chain fructose and tautomers at 303 K
        K_bp303 = 59.2;
        K_bf303 = 26.4;
        K_ap303 = 0.6;
        K_af303 = 6.4;
        
        %Equilibirium constants b/t open-chain fructose and tautomers as a 
        %function of temperature
        K_bp = K_bp303*exp(-(delH_bp/R)*(1/T_K-1/303));
        K_bf = K_bf303*exp(-(delH_bf/R)*(1/T_K-1/303));
        K_ap = K_ap303*exp(-(delH_ap/R)*(1/T_K-1/303));
        K_af = K_af303*exp(-(delH_af/R)*(1/T_K-1/303));
        
        %Furanose fraction at equilibirum as a function of temperature
        phi_f = (K_af+K_bf)/(1+K_af+K_bf+K_ap+K_bp);
        
        %ACTIVATIONS ENERGIES FOR RXN1,RXN2,...,RXN5
        Ea = [127 133 97 64 129]*10^3; %[J/mol]
        
        %NTURAL LOG OF RXN RATE CONSTANTS AT 381 K FOR RXN1,RXN2,...,RXN5
        lnk381 = [1.44 -4.22 -3.25 -5.14 -4.92];
 
        %RXN RATE CONSTANTS FOR RXN1,RXN2,...,RXN5 AS A FUNCTION OF 
        %TEMPERATURE
        k(1) = exp(lnk381(1)-(Ea(1)/R)*(1/T_K-1/Tm)); %[min^-1.M^-1]
        k(2) = exp(lnk381(2)-(Ea(2)/R)*(1/T_K-1/Tm)); %[min^-1.M^-1]
        k(3) = exp(lnk381(3)-(Ea(3)/R)*(1/T_K-1/Tm)); %[min^-1.M^-1]
        k(4) = exp(lnk381(4)-(Ea(4)/R)*(1/T_K-1/Tm)); %[min^-1.M^-1]
        k(5) = exp(lnk381(5)-(Ea(5)/R)*(1/T_K-1/Tm)); %[min^-1.M^-1]
        
        %RXN RATES FOR THE RXN NETWORK OF FRUCTOSE DEHYDRATION
        %Note: C(1) = Normalized Fructose concentration; C(2) = Normalized
        %HMF concentration; C(3) = Normalized LA concentration; 
        %C(4) = Normalized FA concentration;
        Rxn = zeros(1,5); %[mol/L-min]
        Rxn(1) = k(1)*phi_f*C(1)*C_Hplus/C_H2O;
        %[mol/L-min]
        Rxn(2) = k(2)*C(1)*C_Hplus; %[min^-1]
        Rxn(3) = k(3)*C(2)*C_Hplus; %[min^-1]
        Rxn(4) = k(4)*C(2)*C_Hplus; %[min^-1]
        Rxn(5) = k(5)*C(1)*C_Hplus; %[min^-1]
 
        %SPECIES CONSERVATIVE EQUATIONS
        %Notation: rhs = dC/dt
        rhs(1,1) = (-Rxn(1)-Rxn(2)-Rxn(5)); %Fructose
        rhs(2,1) = (Rxn(1)-Rxn(3)-Rxn(4)); %HMF
        rhs(3,1) = Rxn(3); %LA
        rhs(4,1) = (Rxn(3)+Rxn(5)); %FA
  
    end
end

