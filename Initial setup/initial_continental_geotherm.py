import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# material properties
height = 600e3 # total height in m
upper = 20e3 # upper crust thickness in m
lower = 20e3 # lower crust thickness in m
mantle = 80e3 # mantle thickness in m
astheno = 480e3 # asthenosphere thickness in m


# geotherm parameters
t1 = 273 # ground surface temperature in k
t4 = 1579.66667 # lithosphere-asthenosphere boundary temperature in k
t5 = 1793 # model bottom boundary temperature in k

k1 = 2.5 # thermal conductivity of upper crust in W / (m K)
k2 = 2.5 # thermal conductivity of lower crust in W / (m K)
k3 = 3.0 # thermal conductivity of mantle lithosphere in W / (m K)
k4 = 57.15 # thermal conductivity of asthenosphere in W / (m K)

A1 = 1.0e-6
A2 = 0.4e-6
A3 = 0.02e-6
A4 = 0

q1 = 0.055
q2 = q1-A1 * upper
q3 = q2-A2 * lower
q4 = q3-A3 * mantle


# depth range from 0 to 600 km in 600 data points
depth = np.linspace(0, height, 600) 

def calculate_geotherm(_q1, _q2, _q3, _q4, d):
    """
    calculate temperature at a given depth
    _q1, _q2, _q3, _q4 are the top surface heat flux of the layer
    """
    if d >= (height - upper):
        return t1 + (_q1/k1)*(height-d) - (A1*(height-d)**2) / (2*k1)
    elif (height-upper-lower) <= d < (height-upper):
        t2 = t1 + (_q1/k1)*upper - (A1*upper**2) / (2*k1)
        return t2 + (_q2/k2)*(height-d-upper) - (A2*(height-d-upper)**2) / (2*k2)
    elif (height-upper-lower-mantle) <= d < (height-upper-lower):
        t3 = t1 + (_q1/k1)*upper - (A1*upper**2) / (2*k1) + (_q2/k2)*lower - (A2*lower**2) / (2*k2)
        return t3 + (_q3/k3)*(height-d-upper-lower) - (A3*(height-d-upper-lower)**2) / (2*k3)
    else:
        t4 = t1 + (_q1/k1)*upper - (A1*upper**2) / (2*k1) + (_q2/k2)*lower - (A2*lower**2) / (2*k2) + (_q3/k3)*mantle - (A3*mantle**2) / (2*k3)
        return t4 + (_q4/k4)*(height-d-upper-lower-mantle)

temp = np.array([calculate_geotherm(q1, q2, q3, q4, d) for d in depth])


depth_km = (height-depth) / 1e3 # convert depth to kilometers for the plot
data = np.column_stack((depth_km, temp))
np.savetxt('geotherm.csv', data, delimiter=',', header='Depth (km), Temperature (K)', comments='')


# calculate top temperature of each layer for each model
layers_depths = [20e3, 40e3, 120e3, 600e3]
layer_names = ['Upper crust', 'Lower crust', 'Mantle lithosphere', 'Asthenosphere']


print('Temperature for the top surface of each layer:')
for layer_name, depth_m in zip(layer_names, layers_depths):
    temp_ref = calculate_geotherm(q1, q2, q3, q4, height - depth_m)
    print(f'{layer_name}: {temp_ref:.2f} K')


# plot the depth vs temperature graph
fig, ax = plt.subplots(figsize=(6, 6))

ax.plot(temp, depth_km, color='black', label='Reference')

ax.fill_betweenx(depth_km, min(temp), max(temp), where=(depth_km <= 20), color='#222e63', alpha=0.2, label='Upper crust')
ax.fill_betweenx(depth_km, min(temp), max(temp), where=(depth_km > 20) & (depth_km <= 40), color='#6c2f92', alpha=0.2, label='Lower crust')
ax.fill_betweenx(depth_km, min(temp), max(temp), where=(depth_km > 40) & (depth_km <= 120), color='#8d2542', alpha=0.2, label='Mantle lithosphere')
ax.fill_betweenx(depth_km, min(temp), max(temp), where=(depth_km > 120), color='#fac351', alpha=0.2, label='Asthenosphere')

ax.invert_yaxis()
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

minor_locator_x = AutoMinorLocator(4)
ax.xaxis.set_minor_locator(minor_locator_x)
minor_locator_y = AutoMinorLocator(1)
ax.yaxis.set_minor_locator(minor_locator_y)

ax.set_xlim(min(temp), max(temp))

y_min = 0
y_max = 600

ax.set_ylim(y_max, y_min)
ax.set_yticks([0, 20, 40, 120, 600])

ax.set(xlabel='Temperature (K)', ylabel='Depth (km)')

ax.set_title('Initial continental geotherm', fontweight='bold', fontsize=14)
ax.legend()

plt.tight_layout()
plt.savefig('Initial continental geotherm.png', dpi=300, transparent=True)
plt.show()
