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

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - All air has cleared out of the chamber by the average pressure <br />


# Exit Pressure

![Exit Pressure Equation](https://i.ibb.co/tDyTyK4/Screenshot-2022-12-09-134445.png)

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Exhaust is choked in the throat <br />

Numerically solves the expansion ratio formula for exit pressure across chamber pressure during the burn.


# Thrust

![Thrust Equation](https://i.ibb.co/MVBhtS3/Screenshot-2022-12-09-134554.png)

### Assumptions <br />
 - Flow is isentropic <br />
 - Calorically Perfect Gas <br />
 - Alluminum particles are the only non-gas in the exhaust <br />

Calculates the thrust using the rocket thrust equation for chamber pressure and exit pressure during the burn.


# Total Impulse

![Total Impulse Equation](https://i.ibb.co/QYkVXVH/Screenshot-2022-12-09-134634.png)

### Assumptions <br />
 - None <br />

Multiplies the average thrust by the burn time to get total impulse.


# Specific Impulse

![Specific Impulse Equation]()

### Assumptions <br />
 - None <br />
	
Divides total impulse by the mass of the propellent times gravity to get specific impulse.


# Throat Temperature

### Assumptions <br />
 - Flow is isentropic <br />
 - Ideal gas <br />
 - Exhaust is choked in the throat <br />
 - The entire thickness of the throat is the same temperature <br />
 - The chamber is always at the combustion temperture <br />
	
