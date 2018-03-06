
"""Make_systems!!!"""

def make_system_planet(condition,theta): #in degrees
    """
    Using the condition given creates a system with the planet
    at theta with a velocity 90 degress from it allowing it to orbit
    the origin in the counterclockwise direction.
    """
    unpack(condition)

    theta1 = np.deg2rad(theta*UNITS.degree) #transaltes from degress to radians
    x,y = pol2cart(theta1,orbital_radius)   #gets the x and y position given theta and orbital radius
    vx,vy = pol2cart(theta1+.5*pi*UNITS.radian,orbital_speed) #velocity to orbit
    #print(vx,vy)
    init = State(x=x,y=y,vx=vx,vy=vy) #Staaaaate

    ts = linspace(1,duration,ts_f)

    return System(init=init,mass=mass,radius=radius,ts=ts)

def make_system_titan(condition,system,theta): #the system being that of Saturn
    """
    Does the same thing above but in reference to Saturn
    """
    unpack(condition)
    theta1 = np.deg2rad(theta*UNITS.degree)
    x,y = pol2cart(theta1,orbital_radius)

    x += system.init.x #in reference to sat
    y += system.init.y

    vx,vx = pol2cart(theta1+.5*pi*UNITS.radian,orbital_speed)

    init = State(x=x,y=y,vx=vx,vy=vy)

    return System(init=init,mass=mass,radius=radius)

def make_system_combo(condition,theta):

    unpack(condition)

    mars = Condition(orbital_radius = 227.9e9 * m,
              radius = 3396200 * m,
              mass = .64171e24 * kg,
              orbital_speed = 24080 * m/s,
              duration = set_duration,
              ts_f = num_steps
                )

    sys_mars = make_system_planet(mars,theta)

    rocket = Condition(orbital_radius = 100000 *m + mars.radius,
                dry_mass = 1433788 * kg,
                rho = 0.02 * kg/m**3, #needs function to update
                thrust = 22800000 * N,
                fuel_init = 1000000 * kg,
                duration = set_duration,
                ts_f = num_steps
                  )
    sys_rocket = make_system_rocket(rocket,sys_mars,theta)

    mx,my,mvx,mvy = sys_mars.init
    rx,ry,rvx,rvy = sys_rocket.init

    init = State(rx=rx,ry=ry,rvx=rvx,rvy=rvy,mx=mx,my=my,mvx=mvx,mvy=mvy)

    ts = linspace(1,duration,ts_f)

    return System(init=init,m_mass=mars.mass,r_mass=rocket.dry_mass+fuel_init,ts=ts)

def make_system_rocket(condition,system,theta):#system of Mars
    """
    conditon of rocket
    system of mars
    theta of its position relative to mars

    """
    unpack(condition)

    #print(system.init.vx,system.init.vy)
    mvx = system.init.vx
    mvy = system.init.vy

    theta1 = np.deg2rad(theta*UNITS.degree)
    x,y = pol2cart(theta1,orbital_radius)

    xm = system.init.x
    ym = system.init.y

    x += xm #in reference to Mars
    y += ym
    #print(x,y)
    vx,vy = pol2cart(theta1+.5*pi*UNITS.radian,orbital_velocity(Vector(x,y).dist(Vector(xm,ym))))

    vx = vx + mvx
    vy = vy + mvy

    ts = linspace(1,duration,ts_f)

    init = State(x=x,y=y,vx=vx,vy=vy)#,fuel=fuel_init)
    tick = True
    return System(init=init,mass=mass,radius=radius,tick=tick,dry_mass=dry_mass,ts=ts)

#sys_rocket = make_system_rocket(rocket,sys_mars,180)

"""Gravity!!!"""

def calc_mgrav(vr,vm, mass_rocket): #calculate gravity in reference to Mars
    """
    Given the vectors of two objects in space, find the force of gravity acting upon them.
    Pass two vectors with the x and y positons of the rocket and body in question
    G is the gravitational constant
    vr is the position vector of the rocket
    vm is the position vector of mars
    height is just the distance between the center of the two objects
    returns force vector from the rocket towards the body
    """

    height = vr.dist(vm)
    grav = G * mars.mass * mass_rocket / (height)**2
    #print(height,grav/mass_rocket)

    a = Vector(vm.x-vr.x,vm.y-vr.y) #Creates a vector from the rocket to the object
    x,y = pol2cart(a.angle,grav) #the vector has the angle, and with the force of grav we turn them to cartesian

    f_grav = Vector(x,y) #The vector!
    return f_grav

def calc_satgrav(vr,vsat,mass_rocket): #calculate gravity in reference to Saturn
    height = vr.dist(vsat)
    grav = G * saturn.mass * mass_rocket / (height)**2

    a = Vector(vsat.x-vr.x,vsat.y-vr.y)
    x,y = pol2cart(a.angle,grav)

    f_grav = Vector(x,y)
    return f_grav

def calc_tgrav(vr,vt,mass_rocket): #calculate gravity in reference to Titan
    height = vr.dist(vt)
    grav = G * titan.mass * mass_rocket / (height)**2

    a = Vector(vt.x-vr.x,vt.y-vr.y)
    x,y = pol2cart(a.angle,grav)

    f_grav = Vector(x,y)
    return f_grav


def calc_sgrav(vr,mass):#calculate gravity in reference to the Sun
    height = vr.mag #hah
    grav = G * sun.mass * mass / (height)**2

    x,y = pol2cart(vr.angle+pi*UNITS.radians,grav) #added pi to face in the right direction

    f_grav = Vector(x,y)
    return f_grav

def calc_dv(burn_time,system):
    unpack(system)
    """
    Calculates a change in velocity, taking burn time of thrusters as input
    implies a constant thrust, defined above
    refers to STAGE TWO of Spacex Falcon 9, link to source above
    """
    mi = mass_initial
    new_mass = mi + burn_time*s*flow_rate


    h = fuel_mass/new_mass*I
    print(h)


    """
    Calculates change in mass of fuel over time
    references flow rate from source
    """
    #print(saturn.radius)
    dv = g*h*log(mi/new_mass)
    return dv


def calc_radii(dv):
    #dv *= m/s

    vp  =  mars.orbital_speed + dv


    """
    Takes the change in velocity from the previous function
    Adds change in velocity to speed of mars, gives us speed of rocket as it orbits around the sun
    Puts this velocity into a magic sauce equation to calculate the distance
    from the sun to the perigee and apogee, labelled r1 and r2 respectively
    This is the goal but the unit's aren't working, rip

    So I added something to get the right units trusting that your equation works so it should be fine now - Corey
    """
    #print(vp)
    #print(G)
    #print(mars.radius)
    #print(sun.mass)
    #print

    r1 = mars.orbital_radius
    f = (G * sun.mass*r1 - (vp**2)*(r1**2)) / ((vp**2)*(r1)-2*G*sun.mass)
    i=abs(f)
    print(mars.orbital_speed)
    print(vp)
    print(i)

    r2 = 2*i+r1
    return r1, r2

def orbital_velocity(height): #in meters
    """
    passes a height in meters
    returns magnitude (speed) of orbital velocity at that height
    """
    #height *= m
    v = (G*mars.mass/height)**(1/2)
    return v


def slope_func(state,t,system):
    x, y, vx, vy = state
    unpack(system)

    v = Vector(vx,vy)

    a_thrust = accel_thrust(system,v,t)

    a_grav = calc

    a = a_grav + a_thrust

    return v.x, v.y, a.x, a.y

def calc_radii(dv):
    #dv *= m/s

    vp  =  mars.orbital_speed + dv


    """
    Takes the change in velocity from the previous function
    Adds change in velocity to speed of mars, gives us speed of rocket as it orbits around the sun
    Puts this velocity into a magic sauce equation to calculate the distance
    from the sun to the perigee and apogee, labelled r1 and r2 respectively
    This is the goal but the unit's aren't working, rip

    So I added something to get the right units trusting that your equation works so it should be fine now - Corey
    """
    #print(vp)
    #print(G)
    #print(mars.radius)
    #print(sun.mass)
    #print

    r1 = mars.orbital_radius
    f = (G * sun.mass*r1 - (vp**2)*(r1**2)) / ((vp**2)*(r1)-2*G*sun.mass)
    i=abs(f)
    print(mars.orbital_speed)
    print(vp)
    print(i)

    r2 = 2*i+r1
    return r1, r2

def an(sys_rocket,sys_mars):
    rxs = sys_rocket.results.x
    rys = sys_rocket.results.y
    mxs = sys_mars.results.x
    mys = sys_mars.results.y

    animate2d(rxs,rys,mxs,mys,100)

def combo_slope(state,t,system):

    rx,ry,rvx,rvy,mx,my,mvx,mvy = state

    unpack(system)
    #Make some vectors!
    m_pos = Vector(mx,my)
    m_vel = Vector(mvx,mvy)

    r_pos = Vector(rx,ry)
    r_vel = Vector(rvx,rvy)

    #Mars forces...
    m_f_grav = calc_sgrav(m_pos,1)
    m_a_grav = m_f_grav

    #Rocket forces...
    r_f_grav = calc_mgrav(r_pos,m_pos,mass)
    r_a_grav = r_f_grav/mass

    r_f_sgrav = calc_sgrav(r_pos,mass)
    r_a_sgrav = r_f_sgrav/mass

    a = r_a_grav + r_a_sgrav
    #return delta in each
    return rvx,rvy,a.x,a.y,mvx,mvy,m_a_grav.x,m_a_grav.y


def m_slope(state,t,system):

    x,y,vx,vy = state
    unpack(system)

    pos = Vector(x,y)
    vel = Vector(vx,vy)

    f_grav = calc_sgrav(pos,mass)
    a_grav = f_grav / mass

    #a_grav/m**2

    #print(vx)

    return vx,vy,a_grav.x,a_grav.y


def r_slope(state,t,system):

    x,y,vx,vy = state
    unpack(system)
    pos = Vector(x,y)
    vel = Vector(vx,vy)
    #print(pos.mag)

    mpos = Vector(xpos_mars(t),ypos_mars(t))*m

    f_grav = calc_mgrav(pos,mpos,mass)
    a_grav = f_grav/mass

    f_sgrav = calc_sgrav(pos,mass)
    a_sgrav = f_sgrav/mass

    a = a_sgrav + a_grav

    return vx,vy,a.x,a.y

def animate2d(rxs, rys,mxs,mys, speedup=1):
    """Animate the results of a projectile simulation.

    xs: x position as a function of time
    ys: y position as a function of time

    speedup: how much to divide `dt` by
    """
    # get the time intervals between elements
    ts = mxs.index
    dts = np.diff(ts)
    dts = np.append(dts, 0)

    # decorate the plot
    newfig()
    decorate(xlabel='x position (m)',
             ylabel='y position (m)',
             xlim=[rxs.min(), rxs.max()],
             ylim=[rys.min(), rys.max()],
             legend=False)

    # loop through the values

    for rx,ry,mx,my,dt in zip(rxs,rys,mxs,mys,dts):
        plot(0,0,'yo',update=False)
        plot(mx, my, 'ro', update=True)
        plot(rx,ry,'bo',update=True)
        sleep(dt / speedup)

def animate2d_single(sys):
    """Animate the results of a projectile simulation.

    xs: x position as a function of time
    ys: y position as a function of time

    speedup: how much to divide `dt` by
    """
    xs = sys.results.x
    ys = sys.results.y

    # decorate the plot
    newfig()
    decorate(xlabel='x position (m)',
             ylabel='y position (m)',
             xlim=[xs.min(), xs.max()],
             ylim=[ys.min(), ys.max()],
             legend=False)

    # loop through the values

    for x,y in zip(xs,ys):
        plot(0,0,'yo',update=False)
        plot(x,y,'ro',update=True)
        sleep(.01)
