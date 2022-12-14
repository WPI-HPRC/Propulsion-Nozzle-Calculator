
%% Run this after running the main script, it does not work on its own

%% Inputs
T_a = 40; % Temperature of ambient air (F)
K_n = 15; % Thermal conductivity of nozzle material (Wm^-1K^-1)

tube_length = 80; % Length of pressure transducer tube (in)
tube_ID = 0.25; % Inner diameter of pressure transducer tube (in)
emissivity_tube = 0.03; % Emmissivity of tube

%% Conversions
T_a = (5/9)*(T_a-32) + 273.15; % Temperature of ambient air (K)
tube_length = tube_length*0.0254; % Length of pressure transducer tube (m)
tube_ID = tube_ID*0.0254; % Inner diameter of pressure transducer tube (m)

% RK4 Setup
tube_dt = 0.000001;
xRec = zeros(10000000,2);
xCurr = [T_0*(1+((k-1)/2))^-1,0];
xRec(1,:) = xCurr;

% RK4 Loop, uses the distance along the tube to determine when it should
% stop
i = 2;
while(xRec(i-1,2)<tube_length) 

    k1=tube_temperature_dynamics(xCurr,R,k,K_n,T_0,P_avg,M_p,k_b,N_A,d_star,tube_ID/2,emissivity_tube,sigma,T_a)*tube_dt;
    k2=tube_temperature_dynamics(xCurr+1/2*k1,R,k,K_n,T_0,P_avg,M_p,k_b,N_A,d_star,tube_ID/2,emissivity_tube,sigma,T_a)*tube_dt;
    k3=tube_temperature_dynamics(xCurr+1/2*k2,R,k,K_n,T_0,P_avg,M_p,k_b,N_A,d_star,tube_ID/2,emissivity_tube,sigma,T_a)*tube_dt;
    k4=tube_temperature_dynamics(xCurr+k3,R,k,K_n,T_0,P_avg,M_p,k_b,N_A,d_star,tube_ID/2,emissivity_tube,sigma,T_a)*tube_dt;
    xCurr=xCurr+1/6*k1+1/3*k2+1/3*k3+1/6*k4;

    xRec(i,:) = xCurr;
    i = i+1;
end

% Output
T_transducer_final = xRec(i-1,1);

fprintf("Temperature at end of tube: %4.0fK\n\n",T_transducer_final);


function xDot = tube_temperature_dynamics(xCurr,R,k,K_n,T0,P0,M_p,k_b,N_A,d,t_n,emissivity,sigma,T_a)
    T = xCurr(1); 

    Cp = R/(k-1);
    Pr = (2/5)*k;
    mu = Pr*(K_n/Cp);

    v = (k*R*T)^(1/2);
    rho0 = (P0*M_p)/(k_b*T0*N_A);
    rho = rho0*(1+((k-1)/2))^((k-1)/2);
    Re = (rho*v*d)/mu;

    St = 0.023*Re^(-1/5)*Pr^(-0.067);
    h = St*rho*v*Cp;
    q = -((emissivity*sigma*T^4) + (h*(T - T_a)));
    TDot = q/(Cp*rho*t_n);
    pDot = sqrt(k*R*T);

    xDot = [TDot,pDot];
end

