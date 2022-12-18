
%% Run this after running the main script, it does not work on its own

%% Inputs
T_a = 40; % Temperature of ambient air (F)
K_n = 15; % Thermal conductivity of nozzle material (Wm^-1K^-1)
C_n = 500; % Specific heat of nozzle material (m^2s^-2K^-1 | Jkg^-1K^-1)
rho_n = 7500;  % Density of nozzle material (kgm^-3)
t_n = 0.1; % Thickness of nozzle (in);
emissivity = 0.5; % Emissivity of nozzle material

%% Constants
h_a = 10; % Heat transfer coeffient of ambient air (Wm^-2K^-1)

%% Conversions
T_a = (5/9)*(T_a-32) + 273.15; % Temperature of ambient air (K)
t_n = t_n*0.0254; % Thickness of nozzle (m);

% RK4 Setup
xCurr = T_a;

TRec = zeros(length(t),1);
T(1) = xCurr;

% RK4 Loop
i = 1;
while(t(i)<burn_time)

    k1=temperature_dynamics(xCurr,R,kRec(i),K_n,T_0,PRec(i),M_p,k_b,N_A,d_star,C_n,rho_n,t_n,emissivity,sigma,h_a,T_a)*dtRec(i);
    k2=temperature_dynamics(xCurr+1/2*k1,R,kRec(i),K_n,T_0,PRec(i),M_p,k_b,N_A,d_star,C_n,rho_n,t_n,emissivity,sigma,h_a,T_a)*dtRec(i);
    k3=temperature_dynamics(xCurr+1/2*k2,R,kRec(i),K_n,T_0,PRec(i),M_p,k_b,N_A,d_star,C_n,rho_n,t_n,emissivity,sigma,h_a,T_a)*dtRec(i);
    k4=temperature_dynamics(xCurr+k3,R,kRec(i),K_n,T_0,PRec(i),M_p,k_b,N_A,d_star,C_n,rho_n,t_n,emissivity,sigma,h_a,T_a)*dtRec(i);
    xCurr=xCurr+1/6*k1+1/3*k2+1/3*k3+1/6*k4;

    TRec(i) = xCurr;

    i = i+1;
end

% Output
figure()
plot(t(1:i-1),TRec(1:i-1))
title("Temperature of Nozzle over Burn")
xlabel('Time (s)', 'FontSize', 11)
ylabel('Temperature (K)', 'FontSize', 11)

fprintf("Max Temperature of Nozzle: %4.0fK\n\n",max(TRec));


function TDot = temperature_dynamics(TCurr,R,k,K_n,T0,P0,M_p,k_b,N_A,d,C_n,rho_n,t_n,emissivity,sigma,h_a,T_a)
    Cp = R/(k-1);
    Pr = (2/5)*k;
    mu = Pr*(K_n/Cp);

    T = T0*(1+((k-1)/2))^-1;
    v = (k*R*T)^(1/2);
    rho0 = (P0*M_p)/(k_b*T0*N_A);
    rho = rho0*(1+((k-1)/2))^((k-1)/2);
    Re = (rho*v*d)/mu;

    St = 0.023*Re^(-1/5)*Pr^(-0.067);
    h = St*rho*v*Cp;
    QDot = h*(T-TCurr)-((emissivity*sigma*TCurr^4) + (h_a*(TCurr - T_a)));
    TDot = QDot/(C_n*rho_n*t_n);
end

