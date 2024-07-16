import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# material properties
height = 600e3 # total model height in m
upper = 20e3 # upper crust thickness in m
lower = 20e3 # lower crust thickness in m
mantle = 80e3 # mantle lithosphere thickness in m
astheno = 480e3 # asthenosphere thickness in m

upper_rho = 2800 # upper crust density in kg m^(-3)
lower_rho = 2900 # lower crust density in kg m^(-3)
mantle_rho = 3250 # mantle lithosphere density in kg m^(-3)
astheno_rho = 3300 # asthenosphere density in kg m^(-3)

g = 9.81 # gravitational acceleration in m
R = 8.3144626 # ideal gas constant in J mol^(-1) K^(-1)

epsilon = 1.e-15 # reference strain rate in s^(-1)
cohesion = 20e6 # cohesion in Pa
phi = np.radians(30) # internal angle of friction


# dislocation creep parameters
upper_disl_n = 4 # stress exponent of upper crust
upper_disl_A = 8.57e-28 # material constant of upper crust
upper_disl_Q = 223.e3 # activation energy of upper crust
upper_disl_V = 0 # activation volume of upper crust

lower_disl_n = 3 # stress exponent of lower crust
lower_disl_A = 7.13e-18 # material constant of lower crust
lower_disl_Q = 345.e3 # activation energy of lower crust
lower_disl_V = 0 # activation volume of lower crust

mantle_disl_n = 3.5 # stress exponent of mantle lithosphere
mantle_disl_A = 6.52e-16 # material constant of mantle lithosphere
mantle_disl_Q = 530.e3 # activation energy of mantle lithosphere
mantle_disl_V = 18.e-6 # activation volume of mantle lithosphere

astheno_disl_n = 3.5 # stress exponent of asthenosphere
astheno_disl_A = 5.33e-19 # material constant of asthenosphere
astheno_disl_Q = 480.e3 # activation energy of asthenosphere
astheno_disl_V = 11.e-6 # activation volume of asthenosphere


# diffusion creep parameters
grain_size = 5.e-3 # grain size in m

upper_diff_A = 5.97e-19
upper_diff_Q = 223.e3
upper_diff_V = 0
upper_diff_m = 3 # grain size exponent of upper crust

lower_diff_A = 2.99e-25
lower_diff_Q = 159.e3
lower_diff_V = 0
lower_diff_m = 2 # grain size exponent of lower crust

mantle_diff_A = 2.27e-15
mantle_diff_Q = 375.e3
mantle_diff_V = 10.e-6
mantle_diff_m = 3 # grain size exponent of mantle lithosphere

astheno_diff_A = 1.50e-18
astheno_diff_Q = 335.e3
astheno_diff_V = 4.e-6
astheno_diff_m = 3 # grain size exponent of asthenosphere


def disl_power_law(eps, n, A, Q, P, V, T):
    """
    calculate dislocation power law creep
    eps : strain rate in s^(-1)
    n : stress exponent
    A : material constant in P^(-n) s^(-1)
    Q : activation energy [i mol*^(-1)]
    T : temperature [K]
    """

    disl_stress_diff = (eps / A)**(1 / n) * np.exp((Q + P*V) / (n * R * T))
    return disl_stress_diff / 1e6 # convert to MPa


def diff_power_law(eps, d, m, A, Q, P, V, T):
    """
    calculate diffusion power law creep
    eps : strain rate in s^(-1)
    d : grain size
    m : grain size exponent
    A : material constant in Pa^(-1) s^(-1)
    Q : activation energy in J mol^(-1)
    P : pressure [Pa]
    V : activation volume 
    T : temperature [K]
    """

    diff_stress_diff = (eps / A) * np.exp((Q + P*V) / (R * T)) * d**(m)
    return diff_stress_diff / 1e6 # convert to MPa


# load geotherm data from CSV
data = np.genfromtxt('geotherm.csv', delimiter=',', skip_header=1)
depth_km = data[:, 0]
depth = depth_km * 1.e3 # depth in m
temp = data[:, 1] # temperature in K

plastic_diff_stress = np.zeros(len(depth_km)) # differential stress for plastic deformation
disl_diff_stress = np.zeros(len(depth_km)) # differential stress for dislocation creep
diff_diff_stress = np.zeros(len(depth_km)) # differential stress for diffusion creep

density = np.zeros(len(depth_km))
pressure = np.zeros(len(depth_km))


# plot the yield strength envelope
fig, ax = plt.subplots(figsize=(6, 6))

for i, d in enumerate(depth):
    if d <= upper: # upper crust
        density[i] = upper_rho # density variations within layers are small enough to be neglected
        pressure[i] = g * density[i] * d
        plastic_diff_stress[i] = (cohesion * np.cos(phi) + pressure[i] * np.sin(phi)) / 1e6  # in MPa
        disl_diff_stress[i] = disl_power_law(epsilon, upper_disl_n, upper_disl_A, upper_disl_Q, pressure[i], upper_disl_V, temp[i])
        diff_diff_stress[i] = diff_power_law(epsilon, grain_size, upper_diff_m, upper_diff_A, upper_diff_Q, pressure[i], upper_diff_V, temp[i])
    elif d <= (upper + lower): # lower crust
        density[i] = lower_rho
        pressure[i] = g * (density[i] * (d - upper) + upper_rho * upper)
        plastic_diff_stress[i] = (cohesion * np.cos(phi) + pressure[i] * np.sin(phi)) / 1e6
        disl_diff_stress[i] = disl_power_law(epsilon, lower_disl_n, lower_disl_A, lower_disl_Q, pressure[i], lower_disl_V, temp[i])
        diff_diff_stress[i] = diff_power_law(epsilon, grain_size, lower_diff_m, lower_diff_A, lower_diff_Q, pressure[i], lower_diff_V, temp[i])
    elif d <= (upper + lower + mantle): # mantle lithosphere
        density[i] = mantle_rho
        pressure[i] = g * (density[i] * (d - upper - lower) + lower_rho * lower + upper_rho * upper)
        plastic_diff_stress[i] = (cohesion * np.cos(phi) + pressure[i] * np.sin(phi)) / 1e6
        disl_diff_stress[i] = disl_power_law(epsilon, mantle_disl_n, mantle_disl_A, mantle_disl_Q, pressure[i], mantle_disl_V, temp[i])
        diff_diff_stress[i] = diff_power_law(epsilon, grain_size, mantle_diff_m, mantle_diff_A, mantle_diff_Q, pressure[i], mantle_diff_V, temp[i])
    else: # asthenosphere
        density[i] = astheno_rho
        pressure[i] = g * (density[i] * (d - upper - lower - mantle) + mantle_rho * mantle + lower_rho * lower + upper_rho * upper)
        plastic_diff_stress[i] = (cohesion * np.cos(phi) + pressure[i] * np.sin(phi)) / 1e6
        disl_diff_stress[i] = disl_power_law(epsilon, astheno_disl_n, astheno_disl_A, astheno_disl_Q, pressure[i], astheno_disl_V, temp[i])
        diff_diff_stress[i] = diff_power_law(epsilon, grain_size, astheno_diff_m, astheno_diff_A, astheno_diff_Q, pressure[i], astheno_diff_V, temp[i])

min_diff_stress = np.minimum.reduce([plastic_diff_stress, disl_diff_stress, diff_diff_stress])


# plot the yield strength envelope
ax.plot(min_diff_stress, depth_km, color='black', label='Reference')

ax.invert_yaxis()
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

minor_locator_x = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minor_locator_x)
minor_locator_y = AutoMinorLocator(1)
ax.yaxis.set_minor_locator(minor_locator_y)

x_min = 0
x_max = 900
y_min = 0
y_max = 200

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_max, y_min) # y-axis is flipped

ax.set_xlabel('Differential stress (MPa)')
ax.set_ylabel('Depth (km)')

# fill the entire plot with white background first
ax.fill_betweenx(depth_km, x_min, x_max, where=(depth_km >= 0), color='#fff', alpha=1)

# shade the layers
ax.fill_betweenx(depth_km, x_min, x_max, where=(depth_km < 20), color='#222e63', alpha=0.2, label='Upper crust')
ax.fill_betweenx(depth_km, x_min, x_max, where=((depth_km >= 20) & (depth_km < 40)), color='#6c2f92', alpha=0.2, label='Lower crust')
ax.fill_betweenx(depth_km, x_min, x_max, where=((depth_km >= 40) & (depth_km < 120)), color='#8d2542', alpha=0.2, label='Mantle lithosphere')
ax.fill_betweenx(depth_km, x_min, x_max, where=((depth_km > 120)), color='#fac351', alpha=0.2, label='Asthenosphere')

ax.legend(loc='lower right')

ax.set_title('Yield strength envelope', fontweight="bold", fontsize=14)

plt.tight_layout()
plt.savefig('py.Yield strength envelope.png', dpi=300, transparent=True)
plt.show()
