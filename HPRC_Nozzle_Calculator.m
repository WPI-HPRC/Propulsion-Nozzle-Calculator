clear all; close all; clc;

%% Inputs

L = 7; % Length of Grain (in)
d_core = 0.75; % Diameter of the Core (in)
d_star = 0.425; % Diameter of Throat (in)
d_case = 2; % Diameter of Casing (in)
M_p = 0.2367; % Molar of Propellent (kgmol^-1)
T_0 = 2773; % Combustion Temperature (K)
rho_p = 1668.474187; % Density of Solid Propellent (kg/m^3)
m_p = 0.557918615; % Mass of Propellent (kg)
a = 3.51398*10^-5; % Burn Coeffient (ms^-1Pa^1n)
n = 0.327392; % Burn Exponent
k_p = 1.21; % Ratio of Specific Heats of Propellant
ihibited_ends = 0; % Number of Inhibited Ends

%% Conversions

L = L*0.0254; % Length of Grain (m)
d_core = d_core*0.0254; % Diameter of the Core (m)
d_star = d_star*0.0254; % Diameter of Throat (m)
d_case = d_case*0.0254; % Diameter of Casing (m)

%% Constants

n_air = 2.7E25; % Number Density of Air (m^-3)
R_bar = 8.3145; % Universal Gas Constant (kgm^-1s^-2K^-1mol^-1)
N_A = 6.02E23; % Avagadro's Number (mol^-1)
k_b = 1.38E-23; % Boltzman Constant (JK^-1)
P_a = 101325; % Atmospheric Pressure (Pa)
g = 9.81; % Acceleration due to Gravity (ms^-2)
M_a = 0.02897; % Molar Mass of Air (kgmol^-1)
k_a = 1.4; % Ratio of Specific Heats of Air

%% Derived Parameters

R = R_bar/M_p; % Specific Gas Constant (m^2s^-2K^-1)
A_star = pi*(d_star/2)^2; % Area of Throat (m^2)
n_p = rho_p/(M_p*N_A); % Number Density of Solid Propellant (m^-3)


%% Chamebr Pressure

t = zeros(1,10000);

N = n_p*(L*pi*(d_core/2)^2);
xCurr = [N;d_core;L]; 
xRec = zeros(length(xCurr),length(t));
PRec = zeros(1,length(t));
PRec(1) = 101325;

for i=2:length(t)
    dt = 0.00065/(1+exp(-0.0021*(i-2000)));
    k1=chamber_pressure_dynamics(xCurr,L,d_case,M_p,T_0,a,n,k_p,ihibited_ends,R,A_star,N_A,k_b)*dt;
    k2=chamber_pressure_dynamics(xCurr+1/2*k1,L,d_case,M_p,T_0,a,n,k_p,ihibited_ends,R,A_star,N_A,k_b)*dt;
    k3=chamber_pressure_dynamics(xCurr+1/2*k2,L,d_case,M_p,T_0,a,n,k_p,ihibited_ends,R,A_star,N_A,k_b)*dt;
    k4=chamber_pressure_dynamics(xCurr+k3,L,d_case,M_p,T_0,a,n,k_p,ihibited_ends,R,A_star,N_A,k_b)*dt;
    xCurr=xCurr+1/6*k1+1/3*k2+1/3*k3+1/6*k4;
    xRec(:,i)=xCurr;
    t(i)=t(i-1)+dt;

    N = xCurr(1); L_g = xCurr(2); d_c = xCurr(3);
    V = (L-L_g)*pi*(d_case/2)^2 + L_g*pi*(d_c/2)^2;
    P_0 = (N/V)*k_b*T_0;

    PRec(i) = P_0;
end


figure(1)
plot(t,PRec)
title("Chamber Pressure over Burn")
xlabel('time (s)', 'FontSize', 11)
ylabel('Chamber Pressure (Pa)', 'FontSize', 11)


%% Expansion Ratio



%% Exit Pressure



%% Thrust



%% Total Impulse



%% Specific Impulse



%% Functions

function xDot = chamber_pressure_dynamics(x,L,d_case,M_e,T_0,a,n,k_p,ihib,R,A_star,N_A,k_b)
    N = x(1); d_c = x(2); L_g = x(3);
    
    V = (L-L_g)*pi*(d_case/2)^2 + L_g*pi*(d_c/2)^2;
    P_0 = (N/V)*k_b*T_0;
    r = a*P_0^n;
    
    d_cDot = -2*r;
    L_gDot = -(2-ihib)*r;
    NDot = n*change_in_volume(L_g,d_c,r,d_case,ihib)-(N_A/M_e)*mass_flow_rate_out(A_star,P_0,R,T_0,k_p);

    xDot = [NDot;d_cDot;L_gDot];
end

function mDot = mass_flow_rate_out(A_star,P,R,T_0,k)
    P_star = P*((2/(k+1))^(k/(k-1)));
    T_star = T_0*(2/(k+1));

    mDot = A_star*P_star*sqrt(k/(R*T_star))*(((k+1)/2)^((k+1)/(2*(1-k))));
end

function vDot = change_in_volume(Lg,d_core,r,d_case, ihib)
    delta_A_core = pi*(((d_core/2)+r)^2-(d_core/2)^2);
    delta_L = (2-ihib)*r;
 
    vDot = delta_L*pi*((d_case/2)^2-((d_core/2)+r)^2) + Lg*delta_A_core;
end



