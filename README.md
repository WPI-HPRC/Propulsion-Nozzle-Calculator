# Propulsion-Nozzle-Calculator

# Chamber Pressure

Assumptions
	The gasses in the chamber perfectly mix instantaneously
    The entire burn area instantly starts combusting at the start of the burn
    The chamber is always at the combustion temperture
    Constant ratio of specific heats for each species
    All propellent instantly becomes gaseous when its burned
    Flow is isentropic
    Ideal gas
    Exhaust is choked in the throat
    The ends maintain a square corner
    Burn rate is the same everywhere
	
Uses RK4 to numerically solve a system of differential equations for the number of particles in the chamber, core diameter, and length of grain. Solves for chamber pressure from number of particles in the chamber using the ideal gas law.
	
# Expansion Ratio

Assumptions
    Flow is isentropic
    Ideal gas
    All air has cleared out of the chamber by the average pressure

Calculates the expansion ratio based on the average chamber and setting the exit preesure to atmospheric pressure.

# Exit Pressure

Assumptions
    Flow is isentropic
    Ideal gas
    Exhaust is choked in the throat

Numerically solves the expansion ratio formula for exit pressure across chamber pressure during the burn.

# Thrust

Assumptions
    Flow is isentropic
    Calorically Perfect Gas
    Alluminum particles are the only non-gas in the exhaust

Calculates the thrust using the rocket thrust equation for chamber pressure and exit pressure during the burn.

# Total Impulse

Assumptions
    None

Multiplies the average thrust by the burn time to get total impulse.
	
# Specific Impulse

Assumptions
    None
	
Divides total impulse by the mass of the propellent times gravity to get specific impulse.


# Throat Temperature

Assumptions
    Flow is isentropic
    Ideal gas
    Exhaust is choked in the throat
    The entire thickness of the throat is the same temperature
    The chamber is always at the combustion temperture
	
