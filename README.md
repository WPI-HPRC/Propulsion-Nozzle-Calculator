# *Propulsion Nozzle Calculator*

To use the nozzle calculator start by putting in the the input parameters, the exit pressure can be calculated for the optimal value for unchanging atmospheric pressure or can be set to manually be entered, then the chamber preesure and thrust are calculated over the course of the burn as well as the total and specific impulse. <br />
After the main script has been ran, the scripts for casing temperature, throat temperature, or pressure transducer temperature can be ran. Throat temperature should only be ran if a non-ablative material is used(eg. not graphite). The pressure transducer temperature is only important during static fire tests where the pressure is recoded using a tube into the forward closure and should be a large overestimate compared to the real temperature.

# Symbols

$$
\begin{aligned*}
& L: \text{Length of Chamber} \\
& d_{chamber}: \text{Diamter of Chamber} \\
& d_c: \text{Diamter of Core} \\
& L_g: \text{Length of Grain} \\
& e_{inhib}: \text{Number of Inhibited Ends} \\
& A^*: \text{Area of Throat} \\
& M_p: \text{Molar Mass of Propellant} \\
& T_0: \text{Combustion Temperature} \\
& a: \text{Burn Rate Coeffeient} \\
& n: \text{Burn Rate Exponent} \\
& \gamma_p: \text{Ratio of Specific Heats of Propellant} \\
& \gamma_a: \text{Ratio of Specific Heats of Air} \\
& \gamma: \text{Ratio of Specific Heats of Exhaust} \\
& N_p: \text{Number of Particles of Propellant} \\
& N_a: \text{Number of Particles of Air} \\
& N: \text{Total Number of Particles} \\
& N_A: \text{Avagadros Number} \\
& V: \text{Volume of Exhaust Area} \\
& k_b: \text{Boltzmann Constant} \\
& R: \text{Specific Gas Constant} \\
& n_p: \text{Number Density of Propellant} \\
& r: \text{Burn Rate} \\
& P_a: \text{Atmospheric Pressure} \\
& P_e: \text{Exit Pressure} \\
& v_e: \text{Exit Velocity} \\
& F: \text{Thrust Force} \\
& A_e: \text{Exit Area} \\
& I_t: \text{Total Impulse} \\
& F_{avg}: \text{Average Thrust} \\
& t_{burn}: \text{Burn Time} \\
& I_{sp}: \text{Specific Impulse} \\
& m: \text{Mass of Propellant} \\
& g: \text{Standard Gravity} \\
& 
\end{aligned*}
$$

# Chamber Pressure

Uses RK4 to numerically solve a system of differential equations for the number of particles in the chamber, core diameter, and length of grain until the core diamter is larger than the chamber diamter. Solves for chamber pressure from number of particles in the chamber using the ideal gas law.

$$
\begin{aligned*}
& \gamma=\frac{N_p \gamma_p+N_a \gamma_a}{N_p+N_a} \\
& V=\pi L_g \left(\frac{d_c}{2}\right)^2+\pi\left(L-L_g\right)\left(\frac{d_{chamber}}{2}\right)^2 \\
& P_0=\frac{N}{V} \gamma_b T_0 \\
& r=a\left(P_0\right)^n \\
& \dot{V}=(2-e_{inhib})(\pi r)\left[\left(\frac{d_{chamber}}{2}\right)^2-\left(\frac{d_c}{2}+r\right)^2\right]+\pi L_g \left[\left(\frac{d_c}{2}+r\right)^2-\left(\frac{d_c}{2}\right)^2\right] \\
& \dot{m}=A^* P_o \sqrt{\frac{\gamma}{R T_0}}\left(\frac{\gamma+1}{2}\right)^{\frac{\gamma+1}{2(1-\gamma)}} \\
& \underline{\dot{x}}=\left[\begin{array}{c}
\dot{N}_p \\
\dot{N}_a \\
\dot{d}_c \\
\dot{L}_g
\end{array}\right]=\left[\begin{array}{c}
n_p \dot{V}-\frac{N_A}{M_p}\left(\frac{N_p}{N_p+N_a}\right) \dot{m} \\
-\frac{N_A}{M_p}\left(\frac{N_p}{N_p+N_a}\right) \dot{m} \\
2 r \\
(e_{inhib}-2) r
\end{array}\right] \\
&
\end{aligned*}
$$

### Assumptions <br />
 - The gasses in the chamber perfectly mix instantaneously <br />
 - The entire burn area instantly starts combusting at the start of the burn <br />
 - The chamber is always at the combustion temperture <br />
 - Constant ratio of specific heats for each species <br />
 - All propellent instantly becomes gaseous when its burned <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Exhaust is choked in the throat <br />
 - The ends maintain a square corner <br />
 - Burn rate is the same everywhere <br />
	

# Expansion Ratio

Calculates the expansion ratio based on the average chamber and setting the exit preesure to atmospheric pressure. The equation is for the inverse expansion ratio, so we invert it before its used, and we also calulate the exit area and exit diameter from the expansion ratio.

$$
\frac{A^*}{A_e}=\left(\frac{\gamma+1}{2}\right)^{\frac{1}{\gamma-1}}\left(\frac{P_a}{P_0}\right)^{\frac{1}{\gamma}} \sqrt{\left(\frac{\gamma+1}{\gamma-1}\right)\left(1-\left(\frac{P_a}{P_0}\right)^{\frac{\gamma-1}{\gamma}}\right)}
$$

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - All air has cleared out of the chamber by the average pressure <br />


# Exit Pressure

Numerically solves the expansion ratio formula for exit pressure across chamber pressure during the burn. The series converges, so exit pressure is calculated by repeatedly plugging it in until it is close enough to the true exit pressure.

$$
P_e=\left(\frac{A^*}{A_e}\right)^\gamma\left(\frac{\gamma+1}{2}\right)^{\frac{\gamma}{1-\gamma}}\left(\frac{\gamma-1}{\gamma+1}\right)^{\frac{\gamma}{2}} P_0\left(1-\left(\frac{P_e}{P_0}\right)^{\frac{\gamma-1}{\gamma}}\right)^{-\frac{\gamma}{2}}
$$

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Exhaust is choked in the throat <br />


# Thrust

Calculates the thrust using the rocket thrust equation for chamber pressure and exit pressure during the burn. The mass flow rate is shown here, but it is saved during the chamber pressure calculation.

$$
\begin{aligned*}
& \dot{m}=A^* P_o \sqrt{\frac{\gamma}{R T_0}}\left(\frac{\gamma+1}{2}\right)^{\frac{\gamma+1}{2(1-\gamma)}} \\
& v_e=\sqrt{\frac{2 \gamma}{\gamma-1} R T_0\left(1-\left(\frac{P_0}{P_a}\right)^{\frac{1-\gamma}{\gamma}}\right)} \\
& F=\dot{m} v_e+\left(P_e-P_a\right) A_e \\
&
\end{aligned*}
$$

### Assumptions <br />
 - Flow is isentropic <br />
 - Calorically Perfect Gas <br />
 - Alluminum particles are the only non-gas in the exhaust <br />


# Total Impulse

Multiplies the average thrust by the burn time to get total impulse. Average thrust is calculated by numerically integrating thrust over time. Burn time is calculated by taking the final time of the chamber pressure RK4. 

$$
I_t=F_{\text {avg }} t_{\text {burn }}
$$

### Assumptions <br />
 - None <br />


# Specific Impulse

Divides total impulse by the mass of the propellent times gravity to get specific impulse.

$$
I_{s p}=\frac{I_t}{m g}
$$

### Assumptions <br />
 - None <br />
	

# Throat Temperature

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Exhaust is choked in the throat <br />
 - The entire thickness of the throat is the same temperature <br />
 - The chamber is always at the combustion temperture <br />
	
	
# Casing Temperature

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Only one inner liner is used <br />
 - No exhuast gas is able to get around the inner liner  <br />
 - The entire thickness of the phenolic liner is at the same temperature  <br />
 - The entire thickness of the casing is at the same temperature  <br />
 - The chamber is always at the combustion temperture <br />
 
 
 # Pressure Transducer Temperature

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Gas is flowing through the tube <br />
 - The tube is the same temperature as the gas in it <br />
 - Steady state <br />