clear all; close all; clc;

%% Assumptions

% Chamber Pressure
%   The gasses in the chamber perfectly mix instantaneously
%   The entire burn area instantly starts combusting at the start of the burn
%   The chamber is always at the combustion temperture
%   Constant ratio of specific heats for each species
%   All propellent instantly becomes gaseous when its burned
%   Flow is isentropic
%   Ideal gas
%   Exhaust is choked in the throat
%   The ends maintain a square corner
%   Burn rate is the same everywhere

% Expansion Ratio
%   Flow is isentropic
%   Ideal gas
%   All air has cleared out of the chamber by the average pressure

% Exit Pressure
%   Flow is isentropic
%   Ideal gas
%   Exhaust is choked in the throat

% Thrust
%   Flow is isentropic
%   Calorically Perfect Gas
%   Alluminum particles are the only non-gas in the exhaust

% Total Impulse
%   None

% Specific Impulse
%   None

%% Inputs

L = 7; % Length of Casing (in)
L_g = 7; % Length of Grain (in)
d_core = 0.75; % Diameter of the Core (in)
d_star = 0.425; % Diameter of Throat (in)
d_case = 2; % Diameter of Casing (in)
M_p = 0.02367; % Molar of Propellent (kgmol^-1)
T_0 = 2773; % Combustion Temperature (K)
rho_p = 1668.474187; % Density of Solid Propellent (kg/m^3)
a = 3.51398*10^-5; % Burn Coeffient (ms^-1Pa^1n)
n = 0.327392; % Burn Exponent
k_p = 1.21; % Ratio of Specific Heats of Propellant
theta_c = 34.5; % Nozzle Converge Angles (degrees)
theta_d = 17.5; % Nozzle Diverge Angles (degrees)
ihibited_ends = 0; % Number of Inhibited Ends
c_s = 897; % Specific heat of solid particles in exhaust (m^2s^-2K^-1 | Jkg^-1K^-1)
beta = 0.075; % Mass fraction of solid particlesin exhaust

%% Conversions

L = L*0.0254; % Length of Casing (in)
L_g = L_g*0.0254; % Length of Grain (m)
d_core = d_core*0.0254; % Diameter of the Core (m)
d_star = d_star*0.0254; % Diameter of Throat (m)
d_case = d_case*0.0254; % Diameter of Casing (m)
theta_c = theta_c*(pi/180); % Nozzle Converge Angles (radians)
theta_d = theta_d*(pi/180); % Nozzle Diverge Angles (radians)

%% Constants

n_a = 2.7E25; % Number Density of Air (m^-3)
R_bar = 8.3145; % Universal Gas Constant (kgm^-1s^-2K^-1mol^-1 | Jmol^-1K^-1)
N_A = 6.02E23; % Avagadro's Number (mol^-1)
k_b = 1.38E-23; % Boltzman Constant (JK^-1)
P_a = 101325; % Atmospheric Pressure (Pa)
g = 9.81; % Acceleration due to Gravity (ms^-2)
M_a = 0.02897; % Molar Mass of Air (kgmol^-1)
k_a = 1.4; % Ratio of Specific Heats of Air

%% Derived Parameters

R = (1-beta)*(R_bar/M_p); % Specific Gas Constant (m^2s^-2K^-1 | Jkg^-1K^-1)
A_star = pi*(d_star/2)^2; % Area of Throat (m^2)
n_p = (rho_p*N_A)/M_p; % Number Density of Solid Propellant (m^-3)
m_p = rho_p*(L_g*pi*((d_case/2)^2-(d_core/2)^2)); % Mass of Propellent (kg)

%% Settings

tMaxSteps = 100000; % Maximum Amount of Steps for Chamber Pressure Calculation
h = 0.000065; % Chamber Pressure dt Height Parameter
s = 0.0021; % Chamber Pressure dt Shape Parameter
b = 2000;  % Chamber Pressure dt Location Parameter
Accuracy = 0.0001; % Accuracy of Exit Pressure Calculator
First_Guess = 50662; % First Guess of Exit Pressure Calculator (Pa)

%% Chamebr Pressure

t = zeros(1,tMaxSteps);
dtRec = zeros(1,length(t));

V = (L-L_g)*pi*(d_case/2)^2 + L_g*pi*(d_core/2)^2;
N_a = (P_a*V)/(k_b*T_0);
xCurr = [0;N_a;d_core;L_g]; 
PRec = zeros(1,length(t));
PRec(1) = 101325;
mDotRec = zeros(1,length(t));
kRec = zeros(1,length(t));
kRec(1) = 1.4;

i=2;
while(true)
    dt = h/(1+exp(-s*(i-b)));
    dtRec(i-1) = dt;
    
    k1=chamber_pressure_dynamics(xCurr,L,d_case,M_p,T_0,a,n,k_p,ihibited_ends,R,A_star,N_A,k_b,n_p,k_a,M_a,d_star,theta_c,beta,c_s)*dt;
    k2=chamber_pressure_dynamics(xCurr+1/2*k1,L,d_case,M_p,T_0,a,n,k_p,ihibited_ends,R,A_star,N_A,k_b,n_p,k_a,M_a,d_star,theta_c,beta,c_s)*dt;
    k3=chamber_pressure_dynamics(xCurr+1/2*k2,L,d_case,M_p,T_0,a,n,k_p,ihibited_ends,R,A_star,N_A,k_b,n_p,k_a,M_a,d_star,theta_c,beta,c_s)*dt;
    k4=chamber_pressure_dynamics(xCurr+k3,L,d_case,M_p,T_0,a,n,k_p,ihibited_ends,R,A_star,N_A,k_b,n_p,k_a,M_a,d_star,theta_c,beta,c_s)*dt;
    xCurr=xCurr+1/6*k1+1/3*k2+1/3*k3+1/6*k4;
    
    t(i)=t(i-1)+dt;

    N_p = xCurr(1); N_a = xCurr(2); d_c = xCurr(3); L_g = xCurr(4);

    V = (L-L_g)*pi*(d_case/2)^2 + L_g*pi*(d_c/2)^2;

    P_0_a = (N_a/V)*k_b*T_0;
    P_0_p = (N_p/V)*k_b*T_0;
    P_0 = P_0_a + P_0_p;

    PRec(i) = P_0;

    k = (N_a*k_a+N_p*k_p)/(N_a+N_p);
    kRec(i)= two_phase_flow(R,k,beta,c_s);
    mDotRec(i) = mass_flow_rate_out(A_star,P_0,R,T_0,kRec(i));

    if(d_c>d_case&&P_0<P_a)
        break;
    end
    if(i==tMaxSteps)
        fprintf("ERROR: Reached Max Time Steps")
        break;
    end
    i=i+1;
end

figure()
subplot(2,1,1)
plot(t(1:i-1),PRec(1:i-1))
title("Chamber Pressure over Burn")
xlabel('Time (s)', 'FontSize', 11)
ylabel('Chamber Pressure (Pa)', 'FontSize', 11)

iMax = i-1;
burn_time = t(iMax)
P_avg = sum(PRec(1:iMax).*dtRec(1:iMax))/burn_time

%% Expansion Ratio

[expansion_ratio,d_e] = expansion_ratio_calcs(P_avg,P_a,two_phase_flow(R,k_p,beta,c_s),A_star)

A_e = A_star*expansion_ratio;

%% Exit Pressure

P_eRec = zeros(1,length(t));

for u = 1:iMax
    
    P_eRec(u) = exit_pressure_solver(A_star,A_e,kRec(u),PRec(u),First_Guess,Accuracy);

end

P_e_avg = sum(P_eRec(1:iMax).*dtRec(1:iMax))/burn_time

% figure()
% plot(t(1:i-1),P_eRec(1:i-1))
% title("Exit Pressure over Burn")
% xlabel('Time (s)', 'FontSize', 11)
% ylabel('Exit Pressure (Pa)', 'FontSize', 11)

%% Thrust

FRec = zeros(1,length(t));

for u = 1:iMax

    correction = ((1/2)*(1+cos(theta_d)))*0.99*0.995*0.96*0.995;
    FRec(u) = thrust(mDotRec(u),P_eRec(u),P_a,A_e,kRec(u),R,T_0,PRec(u),correction);

end

F_avg = sum(FRec(1:iMax).*dtRec(1:iMax))/burn_time

subplot(2,1,2)
plot(t(1:i-1),FRec(1:i-1))
title("Thrust over Burn")
xlabel('Time (s)', 'FontSize', 11)
ylabel('Thrust (N)', 'FontSize', 11)

%% Total Impulse

I_t = F_avg*burn_time

%% Specific Impulse

Isp = I_t/(m_p*g)

%% Functions

function xDot = chamber_pressure_dynamics(x,L,d_case,M_p,T_0,a,n,k_p,inhib,R,A_star,N_A,k_b,n_p,k_a,M_a,d_star,theta_c,beta,c_s)
    N_p = x(1); N_a = x(2); d_c = x(3); L_g = x(4);

    k = ((N_p*k_p+N_a*k_a)/(N_p+N_a));
    k = two_phase_flow(R,k,beta,c_s);
    
    V = (L-L_g)*pi*(d_case/2)^2 + L_g*pi*(d_c/2)^2;
    P_0 = (N_a/V)*k_b*T_0+(N_p/V)*k_b*T_0;
    r = a*P_0^n;

    vDot = change_in_volume(L_g,d_c,r,d_case,inhib);
    m_pDot = (N_p/(N_p+N_a))*(mass_flow_rate_out(A_star,P_0,R,T_0,k));
    m_aDot = (N_a/(N_p+N_a))*(mass_flow_rate_out(A_star,P_0,R,T_0,k));
    
    if(d_c<d_case)
        N_pDot = n_p*vDot-(N_A/M_p)*m_pDot;
    else
        N_pDot = -(N_A/M_p)*m_pDot;
    end

    N_aDot = -(N_A/M_a)*m_aDot;
    d_cDot = 2*r;
    L_gDot = (inhib-2)*r;    

    xDot = [N_pDot;N_aDot;d_cDot;L_gDot];
end

function mDot = mass_flow_rate_out(A_star,P_0,R,T_0,k)
    mDot = ((A_star*P_0*k)/sqrt(R*T_0))*sqrt((2/(k+1))^((k+1)/(k-1)));
end

function vDot = change_in_volume(L_g,d_core,r,d_case, inhib)
    A_coreDot = pi*(((d_core/2)+r)^2-(d_core/2)^2);
    LDot = (2-inhib)*r;
    A_end = pi*((d_case/2)^2-((d_core/2)+r)^2);
 
    vDot = LDot*A_end + L_g*A_coreDot;
end

function [expan,d_e] = expansion_ratio_calcs(P_0,P_a,k,A_star)
    expan = ((((k+1)/2)^(1/(k-1)))*((P_a/P_0)^(1/k))*sqrt(((k+1)/(k-1))*(1-(P_a/P_0)^((k-1)/k))))^-1;

    A_e = expan*A_star;
    d_e = 2*sqrt(A_e/pi);
end

function P_e = exit_pressure_solver(A_star,A_e,k,P_0,P_e_p,accuracy)
    P_e = ((A_star/A_e)^k)*(((k+1)/2)^(k/(1-k)))*(((k-1)/(k+1))^(k/2))*...
        P_0*(1-((P_e_p/P_0)^((k-1)/k)))^(-k/2);

    if(abs(1-(P_e/P_e_p))>accuracy)
        P_e = exit_pressure_solver(A_star,A_e,k,P_0,P_e,accuracy);
    end
end

function F = thrust(mDot,P_e,P_a,A_e,k,R,T_0,P_0,correction)
    v_e = sqrt(((2*k)/(k-1))*R*T_0*((1-(P_e/P_0)^((k-1)/k))));
    v_e = correction*v_e;
    F = mDot*v_e + A_e*(P_e-P_a);
    if(F<0)
        F = 0;
    end
end

function k = two_phase_flow(R,k,beta,c_s)
    c_p = (k*R)/(k-1);
    c_v = c_p/k;
    
    k = (((1-beta)*c_p+beta*c_s)/((1-beta)*c_v+beta*c_s));
end


