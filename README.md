# **Propulsion Nozzle Calculator**

# Chamber Pressure

Uses RK4 to numerically solve a system of differential equations for the number of particles in the chamber, core diameter, and length of grain. Solves for chamber pressure from number of particles in the chamber using the ideal gas law.



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

Calculates the expansion ratio based on the average chamber and setting the exit preesure to atmospheric pressure.

$\frac{A^*}{A_e}=\left(\frac{\gamma+1}{2}\right)^{\frac{1}{\gamma-1}}\left(\frac{P_a}{P_0}\right)^{\frac{1}{\gamma}} \sqrt{\left(\frac{\gamma+1}{\gamma-1}\right)\left(1-\left(\frac{P_a}{P_0}\right)^{\frac{\gamma-1}{\gamma}}\right)}$

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - All air has cleared out of the chamber by the average pressure <br />


# Exit Pressure

Numerically solves the expansion ratio formula for exit pressure across chamber pressure during the burn.

$P_e=\left(\frac{A^*}{A_e}\right)^\gamma\left(\frac{\gamma+1}{2}\right)^{\frac{\gamma}{1-\gamma}}\left(\frac{\gamma-1}{\gamma+1}\right)^{\frac{\gamma}{2}} P_0\left(1-\left(\frac{P_e}{P_0}\right)^{\frac{\gamma-1}{\gamma}}\right)^{-\frac{\gamma}{2}}$

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Exhaust is choked in the throat <br />


# Thrust

Calculates the thrust using the rocket thrust equation for chamber pressure and exit pressure during the burn.

$F=\dot{m} v_e+\left(P_e-P_a\right) A_e$

### Assumptions <br />
 - Flow is isentropic <br />
 - Calorically Perfect Gas <br />
 - Alluminum particles are the only non-gas in the exhaust <br />


# Total Impulse

Multiplies the average thrust by the burn time to get total impulse.

$I_t=F_{\text {avg }} t_{\text {burn }}$

### Assumptions <br />
 - None <br />


# Specific Impulse

Divides total impulse by the mass of the propellent times gravity to get specific impulse.

$I_{s p}=\frac{I_t}{m g}$

### Assumptions <br />
 - None <br />
	

# Throat Temperature

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Exhaust is choked in the throat <br />
 - The entire thickness of the throat is the same temperature <br />
 - The chamber is always at the combustion temperture <br />
	
