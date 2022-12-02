
%% Run this after running the main script, it does not work on its own

h_l = 0.05;
C_l = 1700;
rho_l = 301.22;
t_l = 0.0013589;
h_c = 152;
C_c = 890;
rho_c = 2710;
t_c = 0.003175;
e_c = 0.1;
h_s = 1360;


xCurr = [T_a,T_a];
xRec = zeros(iMax,2);
xRec(1,:) = xCurr;

for i = 1:iMax

    k1=casing_temperature_dynamics(xCurr,h_l,T_0,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,h_s)*dtRec(i);
    k2=casing_temperature_dynamics(xCurr+1/2*k1,h_l,T_0,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,h_s)*dtRec(i);
    k3=casing_temperature_dynamics(xCurr+1/2*k2,h_l,T_0,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,h_s)*dtRec(i);
    k4=casing_temperature_dynamics(xCurr+k3,h_l,T_0,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,h_s)*dtRec(i);
    xCurr=xCurr+1/6*k1+1/3*k2+1/3*k3+1/6*k4;

    xRec(i,:) = xCurr;
end

T_casing_final = xRec(i,2);

fprintf("Max temperature of casing: %4.0fK\n\n",T_casing_final);


function xDot = casing_temperature_dynamics(xCurr,h_l,T_0,C_l,rho_l,t_l,h_c,h_a,e_c,sigma,C_c,rho_c,t_c,T_a,h_s)
    T_l = xCurr(1); T_c = xCurr(2);

    q = sigma*T_0^4 + h_l*(T_0 - T_l) - h_c*(T_l-T_c);
    TDot_liner = q/(C_l*rho_l*t_l);

    q = h_s*e_c + h_c*(T_l-T_c) - h_a*(T_c - T_a) - e_c*sigma*T_c^4;
    TDot_casing = q/(C_c*rho_c*t_c);

    xDot = [TDot_liner,TDot_casing];
end