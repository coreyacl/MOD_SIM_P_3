

"""First step in any good code:import your libraries!"""
#%matplotlib notebook
from math import *
from modsim import *
from helper import *
#%matplotlib qt5

"""Setup your units"""

m = UNITS.meter
km = UNITS.kilometer
s = UNITS.second
kg = UNITS.kilogram
N = kg * m/(s**2)
kN = kg * m/(s**2) * 1000
G = 6.67408 * (10**-11) * m**3 * kg**-1 * s**-2

"""Setup conditions for bodies"""

set_duration = 200000000
num_steps = 1500
ahh = SweepSeries()

sun = Condition(x = 0 * m, #i'm using conditions instead of systems to reduce the am't of systems we have
             y = 0 * m,
             mass = 1.9891e30 * kg
            )

mars = Condition(orbital_radius = 227.9e9 * m,
              radius = 3396200 * m,
              mass = .64171e24 * kg,
              orbital_speed = 24080 * m/s,
              duration = set_duration,
              ts_f = num_steps
                )

saturn = Condition(orbital_radius = 1.429e12 * m,
                mass = 5.683e26 * kg,
                orbital_speed = 9680 *m/s,
                radius= 58.232e6 * m,
                duration = set_duration,
                ts_f = num_steps
                  )
sad =saturn.orbital_radius #functions down the line not a fan of pulling from function

titan = Condition(radius = 2575e3 *m,
               mass = 1.3455e23 * kg,
               orbital_speed = 5570 *m/s,
               orbital_radius = 1221865e3 * m,
               duration = set_duration,
               ts_f = num_steps
              )

rocket = Condition(orbital_radius = 100000 *m + mars.radius,
                dry_mass = 1420788*kg + 6500 * kg,
                rho = 0.02 * kg/m**3, #needs function to update
                thrust = 22800000 * N,
                fuel_init = 1000000 * kg,
                duration = set_duration,
                ts_f = num_steps
                  )

r_combo = Condition(r_orbital_radius = 100000 *m + mars.radius,
                r_dry_mass = 1433788 * kg,
                m_orbital_radius = 227.9e9 * m,
                m_radius = 3396200 * m,
                m_mass = .64171e24 * kg,
                m_orbital_speed = 24080 * m/s,
                duration = set_duration,
                ts_f = num_steps
)


"""Some other important functions and such..."""

mass_initial = 1420000 * kg
fuel_mass = 1222800 * kg
thrust = System(
thrust_magnitude = 22800461 * N,
flow_rate = -273.6 * kg/s,
mass_initial = 1420000 * kg,
fuel_mass = 1222800 * kg,
burn_time = 375 * s,
mass_final = mass_initial - fuel_mass,
g = 20*(calc_mgrav(Vector(-mars.orbital_radius-(mars.radius+15000*m),0),Vector(-mars.orbital_radius,0*m),1)).mag,#9.81 * m * s**-2
#print(g)
I = 300 * s, #Merlin in vacuum from https://www.reddit.com/r/spacex/comments/3lsm0q/f9ft_vs_f9v11_fuel_mass_flow_rate_isp/
)

a = calc_dv(2401, thrust)
#a

"""Workspace!!!"""

set_duration = 30000000
num_steps = 200

mars = Condition(orbital_radius = 227.9e9 * m,
              radius = 3396200 * m,
              mass = .64171e24 * kg,
              orbital_speed = 24080 * m/s,
              duration = set_duration,
              ts_f = num_steps
                )

rocket = Condition(orbital_radius = 100000 *m + mars.radius,
                dry_mass = 1433788 * kg,
                rho = 0.02 * kg/m**3, #needs function to update
                thrust = 22800000 * N,
                fuel_init = 1000000 * kg,
                duration = set_duration,
                ts_f = num_steps,
                dv = 1#22000*m/s
                  )

r_combo = Condition(duration = set_duration,
                ts_f = num_steps
)



sys_mars = make_system_planet(mars,180)
sys_rocket = make_system_rocket(rocket,sys_mars,180)
run_odeint(sys_mars,m_slope) #run simulation for mars
xpos_mars = interpolate(sys_mars.results.x) #store its position in an interpol
ypos_mars = interpolate(sys_mars.results.y) #one stores the x and the other the y
run_odeint(sys_rocket,r_slope) #tried running the simulation using interpol in the slope function, does not work
                                #yet for some reason. Still need to debug


sys_rocket.results.head
#animate2d_single(sys_rocket)
