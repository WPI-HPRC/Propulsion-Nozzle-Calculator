# Propulsion-Nozzle-Calculator

# Chamber Pressure

Uses RK4 to numerically solve a system of differential equations for the number of particles in the chamber, core diameter, and length of grain. Solves for chamber pressure from number of particles in the chamber using the ideal gas law.

![Chamber Pressure Equation]()

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

![Expansion Ratio Equation](https://i.ibb.co/R3g0bwW/Screenshot-2022-12-09-134326.png)

$\frac{A^*}{A_e}=\left(\frac{k+1}{2}\right)^{\frac{1}{k-1}}\left(\frac{P_a}{P_0}\right)^{\frac{1}{k}} \sqrt{\left(\frac{k+1}{k-1}\right)\left(1-\left(\frac{P_a}{P_0}\right)^{\frac{k-1}{k}}\right)}$

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - All air has cleared out of the chamber by the average pressure <br />


# Exit Pressure

Numerically solves the expansion ratio formula for exit pressure across chamber pressure during the burn.

![Exit Pressure Equation](https://i.ibb.co/tDyTyK4/Screenshot-2022-12-09-134445.png)

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Exhaust is choked in the throat <br />


# Thrust

Calculates the thrust using the rocket thrust equation for chamber pressure and exit pressure during the burn.

![Thrust Equation](https://i.ibb.co/MVBhtS3/Screenshot-2022-12-09-134554.png)

### Assumptions <br />
 - Flow is isentropic <br />
 - Calorically Perfect Gas <br />
 - Alluminum particles are the only non-gas in the exhaust <br />


# Total Impulse

Multiplies the average thrust by the burn time to get total impulse.

![Total Impulse Equation](https://i.ibb.co/QYkVXVH/Screenshot-2022-12-09-134634.png)

### Assumptions <br />
 - None <br />


# Specific Impulse

Divides total impulse by the mass of the propellent times gravity to get specific impulse.

![Specific Impulse Equation]()

### Assumptions <br />
 - None <br />
	

# Throat Temperature

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Exhaust is choked in the throat <br />
 - The entire thickness of the throat is the same temperature <br />
 - The chamber is always at the combustion temperture <br />
	
