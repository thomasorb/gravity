import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches

import mars_atm_model


g0 = 9.807 # m/s²

def compute_mdot(thrust_vac, Ispvac):
    # thrust_vac in N
    # Ispvac in s
    return thrust_vac / (Ispvac * g0) 

# https://en.wikipedia.org/wiki/Ariane_5
# https://core.ac.uk/download/pdf/77231853.pdf
# https://www.ariane.group/wp-content/uploads/2020/06/VULCAIN2.1_2020_04_PS_EN_Web.pdf
# http://www.astronautix.com/
ariane_stage1_EAP_P241 = { # equipped with vulcain 2 engine
    'diameter':5.4, # used to compute drag
    'nozzle_exit_diameter':2.1,
    'Ispvac':432, 
    'dry_mass':14.7e3, 
    'thrust_vac':1390e3
    #'mdot'=323
}

ariane_stage2_ESC_A = { # 
    'diameter':0, # used to compute drag, but already considered in stage 1
    'nozzle_exit_diameter':0.99,
    'Ispvac':446, 
    'dry_mass':4.54e3, 
    'thrust_vac':67e3
    #'mdot'=323
}

ariane_booster_EPC_H173 = {
    'diameter':3.06, # used to compute drag
    'nozzle_exit_diameter':3.06,
    'Ispvac':275, 
    'dry_mass':33e3,
    'thrust_vac':7080e3
}

ariane_double_booster_EPC_H173 = {
    'diameter':3.06, # used to compute drag
    'nozzle_exit_diameter':3.06*2,
    'Ispvac':275, 
    'dry_mass':33e3*2,
    'thrust_vac':7080e3*2
}

ariane_stage1 = ariane_stage1_EAP_P241
ariane_stage2 = ariane_stage2_ESC_A
ariane_booster = ariane_booster_EPC_H173
ariane_double_booster = ariane_double_booster_EPC_H173

class Engine():
    
    def __init__(self, propellant_mass, model=ariane_stage1, stage_number=1):
        
        # default from Vulcain 2: https://en.wikipedia.org/wiki/Vulcain_(rocket_engine) 
        # and https://arc.aiaa.org/doi/10.2514/1.A33363

        self.stage_number = stage_number
        self.started = False
        self.Ae = model['nozzle_exit_diameter']**2 * np.pi / 4
        self.Ispvac = model['Ispvac'] # s
        self.mdot = compute_mdot(model['thrust_vac'], model['Ispvac']) # kg/s
        self.propellant_mass = propellant_mass
        self.dry_mass = model['dry_mass']
    
    def is_started(self): 
        return bool(self.started)
    
    def start(self):
        if self.is_started(): print('warning: engine already started')
        else: self.started = True

    def stop(self):
        if not self.is_started(): print('warning: engine not started')
        else: self.started = False

    def thrust(self, dt, p0=101325):
        if not self.is_started(): return 0
        if self.is_empty(): return 0
        # dt = thrusting time in s
        # p0 = external ambiant pressure in Pa
        
        # get only a fraction of thrust from what's left of propellant
        needed_mass = self.mdot * dt
        ratio = min(1, self.propellant_mass / needed_mass)
        self.propellant_mass -= needed_mass * ratio
        
        # https://en.wikipedia.org/wiki/Rocket_engine_nozzle
        F = self.Ispvac * g0 * self.mdot - self.Ae * p0 # F in N
        
        return F * ratio
    
    def get_mass(self):
        return self.dry_mass + self.propellant_mass
    
    def is_empty(self):
        if self.propellant_mass > 0: return False
        return True


class AtmosphericLayer(object):
    
    def __init__(self, altitude_range, T_params, P_params, R_params, model='exponential', ):


        # h in meters
        
        self.altitude_range = np.array(altitude_range)
        
        if model == 'exponential':
            self.T = lambda h: T_params[0] + T_params[1] * h # degrees Celsius
            self.P = lambda h: P_params[0] + np.exp(P_params[1] * h) # kPa
            self.R = lambda h: np.squeeze(self.P(h) / [.1921 * (self.T(h) + 273.15)]) # kg/m^3
            
        elif model == 'power':
            pass

        elif model == 'function':
            self.T = T_params
            self.P = P_params
            self.R = R_params
        else: raise Exception('bad atmospheric model type, must be exponential or power')

    def plot(self):

        h = np.linspace(*self.altitude_range, 1000)

        fig, axes = plt.subplots(1, 3, figsize=(15,5))
        axes[0].plot(self.R(h), h)
        axes[0].set_title('Density')
        axes[0].set_xlabel(r'$\rho$ (kg/m$^3$)')
        axes[0].set_ylabel('altitude (m)')
        axes[0].grid()

        axes[1].plot(self.T(h), h)
        axes[1].set_title('Temperature')
        axes[1].set_xlabel(r'T ($^\circ$C)')
        axes[1].set_ylabel('altitude (m)')
        axes[1].grid()


        axes[2].plot(self.P(h), h)
        axes[2].set_title('Pressure')
        axes[2].set_xlabel(r'P (kPa)')
        axes[2].set_ylabel('altitude (m)')
        axes[2].grid()

class Atmosphere(object):

    
    def __init__(self):
        self.layers = list()
        self.h_min = None
        self.h_max = None
        self.h_limits = list()

    def add_layer(self, layer):

        self.layers.append(layer)
        self.update()

    def update(self):

        h = list()
        T = list()
        P = list()
        R = list()
        
        for layer in self.layers:
            h_min = np.min(layer.altitude_range)
            h_max = np.max(layer.altitude_range)
            self.h_limits.append(h_max)
            
            if self.h_min is not None:
                self.h_min = min(self.h_min, h_min)
            else:
                self.h_min = h_min
                
            if self.h_max is not None:
                self.h_max = max(self.h_max, h_max)
            else:
                self.h_max = h_max
                

    def get_parameter(self, h, param):

        h = np.atleast_1d(h)
        p = np.full_like(h, np.nan, dtype=float)
        h[h < self.h_min] = self.h_min
        h[h > self.h_max] = self.h_max

        for layer in self.layers:

            ok = (h >= np.min(layer.altitude_range)) * (h <= np.max(layer.altitude_range))
            
            if param == 'T':
                f = layer.T
            elif param == 'P':
                f = layer.P
            elif param == 'R':
                f = layer.R
            else: raise Exception('param must be T, P or R')
                
            p[ok] = f(h[ok])

        if np.any(np.isnan(p)):
            raise Exception('The atmosphere is ill defined, some ranges of altitudes between the layers may not be defined')

        return p


    def get_R(self, h):
        return self.get_parameter(h, 'R')
    def get_T(self, h):
        return self.get_parameter(h, 'T')
    def get_P(self, h):
        return self.get_parameter(h, 'P')
            

    def plot(self):

        def plot_limits(ax):
            for ih in self.h_limits:
                ax.axhline(ih, c='gray', ls='--')
                
        h = np.linspace(self.h_min, self.h_max, 1000)
        
        fig, axes = plt.subplots(1, 3, figsize=(15,5))
        axes[0].plot(self.get_R(h), h)
        axes[0].set_title('Density')
        axes[0].set_xlabel(r'$\rho$ (kg/m$^3$)')
        axes[0].set_ylabel('altitude (m)')
        plot_limits(axes[0])
        axes[0].set_xscale('log')
        axes[0].grid()

        axes[1].plot(self.get_T(h), h)
        axes[1].set_title('Temperature')
        axes[1].set_xlabel(r'T ($^\circ$C)')
        axes[1].set_ylabel('altitude (m)')
        #axes[1].set_xscale('log')
        axes[1].grid()


        axes[2].plot(self.get_P(h), h)
        axes[2].set_title('Pressure')
        axes[2].set_xlabel(r'P (kPa)')
        axes[2].set_ylabel('altitude (m)')
        axes[2].set_xscale('log')
        plot_limits(axes[2])
        axes[2].grid()



class Planet(object):

    
    def __init__(self, mass, radius, atm_model):
        self.mass = mass # kg
        self.radius = radius # m
        self.atm_model = atm_model
        self.G = 6.6743e-11

    def get_g(self, h):
        # return gravitational acceleration in m/s^2
        
        return self.G * self.mass / (self.radius + h)**2

    def get_R(self, h):
        return self.atm_model.get_R(h)

    def get_P(self, h):
        return self.atm_model.get_P(h)

    def plot(self):
        self.atm_model.plot()

        h = np.linspace(self.atm_model.h_min, self.atm_model.h_max, 1000)
        plt.figure()
        plt.plot(self.get_g(h), h)
        plt.grid()
        plt.xlabel('g (m/s$^2$)')
        plt.ylabel('altitude (m)')
        plt.title('Gravitational Acceleration')

    def compute_htheta(self, xy):
        h = np.sqrt(np.sum(xy**2)) - self.radius
        theta = np.arctan2(xy[1], xy[0])
        return h, theta

    def get_PRgh(self, xy, v):
        # xy and v in cartesian coordinates
        
        h, theta = self.compute_htheta(xy)
        g = -self.get_g(h) * np.array([np.cos(theta), np.sin(theta)])
        R = self.get_R(h)
        P = self.get_P(h)
        return P, R, g, h

    def get_circular_orbit_velocity(self, h):

        return np.sqrt(self.G * self.mass / (h + self.radius))

class SpaceCraft(object):

    def __init__(self, mass, p0, v0, a0, area, F0=(0,0), gridsize_min=3, min_dt=0.1, max_dt=1000, engines=None):

        self.mass = mass # payload in kg
        
        self.engines = list()
        if engines is not None:
            for istage, (iengine, ipropellant_mass) in enumerate(engines):
                self.engines.append(Engine(ipropellant_mass, model=iengine, stage_number=istage+1))
            
        self.max_microthrust = 0 # in N
        self.microthrust = 0
        self.F = np.array([0,0]).astype(float) # Force applied (Fx in N, Fy in N)
        self.p = np.array(p0).astype(float) # (x0 in m, y0 in m)
        self.v = np.array(v0).astype(float) # (v0x (m/s), v0y (m/s)
        self.a = np.array(a0).astype(float) # (a0x (m/s^2), a0y (m/s^2))
        self.parachute_tension = np.array([0,0]).astype(float) # (Tx in N, Ty in N)
        self.probe_area = float(area)
        self.area = float(area) # in m^2
        self.gridsize_min = gridsize_min # minimum grid size in m
        self.max_dt = max_dt # max dt in s
        self.min_dt = min_dt # max dt in s

        self.all_F = [F0, ]
        self.all_p = [p0, ]
        self.all_v = [v0, ]
        self.all_a = [a0, ]
        self.all_t = [0, ]
        self.t = 0
        self.all_parachute_tension = [self.parachute_tension, ]
        self.parachute_deployed = False
        self.all_parachute_deployed = [self.parachute_deployed, ]
        self.microthruster_on = False
        self.all_microthruster_on = [self.microthruster_on, ]
        self.all_microthrust = [self.microthrust, ]
        self.thruster_on = False
        self.all_thruster_on = [self.thruster_on, ]
        
    def update(self, P, R, g, Fext=np.array([0.,0.]), dt=None):

        # compute a good dt given the minimum grid size and the maximum dt
        if dt is None:
            V = np.sqrt(np.sum(self.v**2))
            A = np.sqrt(np.sum(self.a**2))
            if A > 0:
                delta = V**2 + 2 * A * self.gridsize_min
                dt = max((-V - np.sqrt(delta)) / A, (-V + np.sqrt(delta)) / A)
            else:
                dt = self.gridsize_min / V
            dt = max(min(dt, self.max_dt), self.min_dt)

        if self.microthruster_on:
            self.microthrust = self.along_v(float(self.max_microthrust))
        else:
            self.microthrust = 0
            self.max_microthrust = 0
            
        # add ext forces
        self.F = np.array(Fext)

        # add drag force
        v = np.sqrt(np.sum(self.v**2))
        drag = 0.5 * R * self.area * v**2
        drag = -self.along_v(drag)
        self.F += drag

        # add g and microthrust forces
        self.F += (self.get_mass()) * g
        self.F += self.microthrust

        # detach empty engines
        for i in range(len(self.engines)):
            if self.engines[i].is_empty():
                print(f'engine {self.engines[0].stage_number} detached')
                del self.engines[i]
                break

        # add engine force
        if len(self.engines) > 0:
            thrust = self.along_v(self.engines[0].thrust(dt, p0=P))
            self.F += thrust
        
        self.a = self.F / (self.get_mass())
        self.v += self.a * dt
        self.p += self.v * dt

        self.all_p.append(np.array(self.p))
        self.all_v.append(np.array(self.v))
        self.all_a.append(np.array(self.a))
        self.all_t.append(float(dt))

        if self.parachute_deployed:
            self.parachute_tension = self.get_mass() * (self.a + g)
        else:
            self.parachute_tension = np.array([0,0]).astype(float)
        self.all_parachute_tension.append(np.array(self.parachute_tension))
        self.all_parachute_deployed.append(np.array(self.parachute_deployed))
            
        self.all_microthruster_on.append(bool(self.microthruster_on))
        self.all_microthrust.append(float(self.max_microthrust))
        self.t += dt
        return dt

    def get_mass(self):
        mass = float(self.mass)
        for iengine in self.engines:
            mass += iengine.get_mass()
        
        return mass
    
    def deploy_parachute(self, area):
        self.area = float(area)
        self.parachute_deployed = True

    def detach_parachute(self):
        if not self.parachute_deployed:
            raise Exception('parachute was not deployed')
        self.area = float(self.probe_area)
        self.parachute_deployed = False

    def start_microthruster(self, max_microthrust):
        self.max_microthrust = max_microthrust
        self.microthruster_on = True

    def stop_microthruster(self):
        self.microthruster_on = False

    def start_thruster(self):
        if len(self.engines) > 0:
            if not self.engines[0].is_started():
                self.engines[0].start()
                print(f'engine {self.engines[0].stage_number} started')
        

    def stop_thruster(self):
        if len(self.engines) == 0:
            raise Exception('No engine attached')
        self.engines[0].stop()

    def along_v(self, a):
        # project a scalar along v vector
        v_theta = np.arctan2(self.v[1], self.v[0])
        return a * np.array([np.cos(v_theta), np.sin(v_theta)])

    def plot_trajectory(self, planet, theta_target=None, overplot=False, altitude_target=None, alpha=0.01, color='tab:blue'):

        if not overplot:
            fig = plt.figure(figsize=(10,10))
            
        xy = np.array(self.all_p)
        colors = np.where(self.all_parachute_deployed, 'tab:green', color)
        colors = np.where(np.array(self.all_microthrust) > 0, 'purple', colors)
        colors = np.where(np.array(self.all_microthrust) < 0, 'red', colors)
        
        plt.scatter(xy[:,0], xy[:,1], marker='.', c=colors, alpha=alpha)

        if not overplot:
            for ih in planet.atm_model.h_limits:
                #cc = matplotlib.patches.Circle((0., 0.), planet.radius + ih, color='gray', fill=False, alpha=0.3)
                #fig.gca().add_patch(cc)
                
                cc = matplotlib.patches.Circle((0., 0.), planet.radius + ih, color='gray', fill=True, alpha=0.1)
                fig.gca().add_patch(cc)
                
            cc = matplotlib.patches.Circle((0., 0.), planet.radius, color='moccasin', fill=True)
            fig.gca().add_patch(cc)

        if theta_target is not None:
            thetas = np.linspace(*theta_target, 100)
            target_xy = [planet.radius * np.cos(np.deg2rad(90-thetas)),
                         planet.radius * np.sin(np.deg2rad(90-thetas))]
            plt.plot(*target_xy, c='red')

        if altitude_target is not None:
            cc = matplotlib.patches.Circle((0., 0.), planet.radius + altitude_target, color='red', fill=False, alpha=0.5)
            fig.gca().add_patch(cc)
                
        
        xlim = np.min(xy[:,0]), np.max(xy[:,0])
        ylim = np.min(xy[:,1]), np.max(xy[:,1])

        if theta_target is not None:
            xlim = min(xlim[0], np.min(target_xy[0])), max(xlim[1], np.max(target_xy[0]))
            ylim = min(ylim[0], np.min(target_xy[1])), max(ylim[1], np.max(target_xy[1]))
        
            
        dx = max((xlim[1] - xlim[0]), (ylim[1] - ylim[0])) * 0.1
        dy = dx


        t = np.cumsum(self.all_t)
        times = np.arange(0, np.size(t), np.size(t)//11)
        for itime in times:
            plt.scatter(xy[itime,0], xy[itime,1], marker='.', alpha=1, color='red')
            plt.text(xy[itime,0] + dx*0.1, xy[itime,1]+dy*0.1, '{:.0f}s'.format(t[itime]))

        plt.xlim(xlim[0] - dx, xlim[1] + dx)
        plt.ylim(ylim[0] - dy, ylim[1] + dy)
            
        if not overplot:
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')
            plt.title('Probe Trajectory')
    
            fig.gca().set_aspect('equal', 'box')
            plt.grid()


    def plot_velocity(self):

        plt.figure()
        plt.plot(np.cumsum(self.all_t), np.sqrt(np.sum(np.array(self.all_v)**2, axis=1)))
        #plt.plot(np.cumsum(self.all_t), np.array(self.all_v))
        #plt.yscale('log')
        plt.xlabel('t (s)')
        plt.ylabel('v (m/s)')
        plt.title('Probe Velocity')
        plt.grid()

    def plot_acceleration(self):

        plt.figure()
        plt.plot(np.cumsum(self.all_t), np.sqrt(np.sum(np.array(self.all_a)**2, axis=1)))
        #plt.yscale('log')
        plt.xlabel('t (s)')
        plt.ylabel(r'a (m/s$^2$)')
        plt.title('Probe Acceleration')
        plt.grid()

    def plot_parachute_tension(self, max_tension=5e3):
        plt.figure()
        #plt.plot(np.cumsum(self.all_t), self.all_parachute_tension)
        plt.plot(np.cumsum(self.all_t), np.sqrt(np.sum(np.array(self.all_parachute_tension)**2, axis=1)))
        plt.axhline(max_tension, c='red')
        #plt.yscale('log')
        plt.xlabel('t (s)')
        plt.ylabel(r'T (N)')
        plt.title('Parachute Tension')
        plt.grid()
        

mars_atm = Atmosphere()
# https://www.grc.nasa.gov/www/k-12/airplane/atmosmrm.html
#mars_atm.add_layer(AtmosphericLayer([0, 7000], [-31, - 0.000998], [.699, -0.00009], [.1921]))
#mars_atm.add_layer(AtmosphericLayer([7000, 20000], [-23.4, - 0.00222 ], [.699, -0.00009], [.1921]))
mars_atm.add_layer(AtmosphericLayer([1, 120000], mars_atm_model.ftemp, mars_atm_model.fpres, mars_atm_model.fdens, model='function'))
Mars = Planet(6.417e23, 3389.5e3, mars_atm)