
%% Run this after running the main script, it does not work on its own

%% Inputs
T_a = 40; % Temperature of ambient air (F)

h_l = 0.05; % Heat transfer coeffient of liner (Wm^-2K^-1)
C_l = 1700; % Specific heat of liner (Jkg^-1K^-1)
rho_l = 301.22; % Density of liner (kgm^-3)
t_l = 0.0535; % Thickness of liner (in)
e_l = 0.95; % Emmisivity of liner

h_c = 152; % Heat transfer coeffient of casing (Wm^-2K^-1)
C_c = 890; % Specific heat of casing (Jkg^-1K^-1)
rho_c = 2710; % Density of casing (kgm^-3)
t_c = 0.125; % Thickness of casing (in)
e_c = 0.1; % Emmisivity of casing

%% Constants
q_s = 1360; % Heat from sun (J/m^2)
h_a = 10; % Heat transfer coeffient of ambient air (Wm^-2K^-1)
t_l = t_l*0.0254; % Thickness of liner (m)
t_c = t_c*0.0254; % Thickness of casing (m)

%% Conversions
T_a = (5/9)*(T_a-32) + 273.15; % Temperature of ambient air (K)

% RK4 Setup
T = T_0;
xCurr = [T_a,T_a];
xRec = zeros(iMax*10,2);
xRec(1,:) = xCurr;

% RK4 Loop
for i = 1:iMax*10
    if(i>iMax) % Checks whether the bunr has ended, when it swaps over the ambient temp inside the chamber
        dt = dtRec(iMax);
        T = T_a;
    else
        dt = dtRec(i);
    end
    

    k1=casing_temperature_dynamics(xCurr,h_l,T,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,q_s,e_l)*dt;
    k2=casing_temperature_dynamics(xCurr+1/2*k1,h_l,T,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,q_s,e_l)*dt;
    k3=casing_temperature_dynamics(xCurr+1/2*k2,h_l,T,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,q_s,e_l)*dt;
    k4=casing_temperature_dynamics(xCurr+k3,h_l,T,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,q_s,e_l)*dt;
    xCurr=xCurr+1/6*k1+1/3*k2+1/3*k3+1/6*k4;

    xRec(i,:) = xCurr;
end

% Output
T_casing_final = xRec(iMax,2);
T_casing_max = max(xRec(:,2));

fprintf("Max temperature of casing during burn: %4.0fK\n",T_casing_final);
fprintf("Max temperature of casing: %4.0fK\n",T_casing_max);


function xDot = casing_temperature_dynamics(xCurr,h_l,T_0,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,q_s,e_l)
    T_l = xCurr(1); T_c = xCurr(2);

    q = e_l*sigma*(T_0^4-T_l^4) + h_l*(T_0 - T_l) - h_c*(T_l-T_c);
    TDot_liner = q/(C_l*rho_l*t_l);

    q = q_s*e_c + h_c*(T_l-T_c) - h_a*(T_c - T_a) - e_c*sigma*(T_c^4-T_a^4);
    TDot_casing = q/(C_c*rho_c*t_c);

    xDot = [TDot_liner,TDot_casing];
end
