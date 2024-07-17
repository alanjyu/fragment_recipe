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

k1 = 2.5 # thermal conductivity of upper crust in W m^(-1) K^(-1)
k2 = 2.5 # thermal conductivity of lower crust in W m^(-1) K^(-1)
k3 = 3.0 # thermal conductivity of mantle lithosphere in W m^(-1) K^(-1)
k4 = 57.15 # thermal conductivity of asthenosphere in W m^(-1) K^(-1)

H1 = 1.0e-6 # rate of heat generation per volume of upper crust in W m^(-3)
H2 = 0.4e-6 # rate of heat generation per volume of lower crust in W m^(-3)
H3 = 0.02e-6 # rate of heat generation per volume of mantle lithosphere in W m^(-3)
H4 = 0 # rate of heat generation per volume of asthenosphere in W m^(-3)

q1 = 0.055 # top surface heat flux of upper crust in W m^(-2)
q2 = q1-H1 * upper # top surface heat flux of lower crust in W m^(-2)
q3 = q2-H2 * lower # top surface heat flux of mantle lithosphere in W m^(-2)
q4 = q3-H3 * mantle # top surface heat flux of asthenosphere in W m^(-2)


# depth range from 0 to 600 km in 600 data points
depth = np.linspace(0, height, 600) 

def calculate_geotherm(_q1, _q2, _q3, _q4, d):
    """
    calculate temperature at a given depth
    _q1, _q2, _q3, _q4 are the top surface heat flux of the layer
    """
    if d >= (height - upper):
        return t1 + (_q1/k1)*(height-d) - (H1*(height-d)**2) / (2*k1)
    elif (height-upper-lower) <= d < (height-upper):
        t2 = t1 + (_q1/k1)*upper - (H1*upper**2) / (2*k1)
        return t2 + (_q2/k2)*(height-d-upper) - (H2*(height-d-upper)**2) / (2*k2)
    elif (height-upper-lower-mantle) <= d < (height-upper-lower):
        t3 = t1 + (_q1/k1)*upper - (H1*upper**2) / (2*k1) + (_q2/k2)*lower - (H2*lower**2) / (2*k2)
        return t3 + (_q3/k3)*(height-d-upper-lower) - (H3*(height-d-upper-lower)**2) / (2*k3)
    else:
        t4 = t1 + (_q1/k1)*upper - (H1*upper**2) / (2*k1) + (_q2/k2)*lower - (H2*lower**2) / (2*k2) + (_q3/k3)*mantle - (H3*mantle**2) / (2*k3)
        return t4 + (_q4/k4)*(height-d-upper-lower-mantle)

temp = np.array([calculate_geotherm(q1, q2, q3, q4, d) for d in depth])


depth_km = (height-depth) / 1e3 # convert depth to kilometers for the plot
data = np.column_stack((depth_km, temp))
np.savetxt('geotherm.csv', data, delimiter=',', header='Depth (km), Temperature (K)', comments='')


# calculate the bottom temperature for each layer
layers_depths = [upper, upper+lower, upper+lower+mantle, height]
layers = ['Upper crust', 'Lower crust', 'Mantle lithosphere', 'Asthenosphere']
bottom_temp = np.zeros(len(layers))

print('Temperature for the bottom surface of each layer:')
for i, depth_m in enumerate(layers_depths):
    bottom_temp[i] = calculate_geotherm(q1, q2, q3, q4, height - depth_m)
    print(f'{layers[i]}: {bottom_temp[i]:.2f} K')


# plot depth versus temperature graph
fig, ax = plt.subplots(figsize=(6, 6))

ax.plot(temp, depth_km, color='black', label='Reference')

ax.invert_yaxis()
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

minor_locator_x = AutoMinorLocator(4)
ax.xaxis.set_minor_locator(minor_locator_x)
minor_locator_y = AutoMinorLocator(1)
ax.yaxis.set_minor_locator(minor_locator_y)

x_min = min(temp)
x_max = max(temp)
y_min = 0
y_max = 600

ax.set_xlim(min(temp), max(temp))
ax.set_xticks([min(temp), bottom_temp[0], bottom_temp[1], bottom_temp[2], bottom_temp[3]])

ax.set_ylim(y_max, y_min)
ax.set_yticks([y_min, upper/1e3, (upper+lower)/1e3, (upper+lower+mantle)/1e3, height/1e3])

ax.set(xlabel='Temperature (K)', ylabel='Depth (km)')

ax.set_title('Initial continental geotherm', fontweight='bold', fontsize=14)

# fill the entire plot with white background first
ax.fill_betweenx(depth_km, min(temp), max(temp), where=(depth_km >= 0), color='#fff', alpha=1)

# shade the layers
ax.fill_betweenx(depth_km, min(temp), max(temp), where=(depth_km <= upper/1e3), color='#222e63', alpha=0.2, label=layers[0])
ax.fill_betweenx(depth_km, min(temp), max(temp), where=(depth_km > upper/1e3) & (depth_km <= (upper+lower)/1e3), color='#6c2f92', alpha=0.2, label=layers[1])
ax.fill_betweenx(depth_km, min(temp), max(temp), where=(depth_km > (upper+lower)/1e3) & (depth_km <= (upper+lower+mantle)/1e3), color='#8d2542', alpha=0.2, label=layers[2])
ax.fill_betweenx(depth_km, min(temp), max(temp), where=(depth_km > (upper+lower+mantle)/1e3), color='#fac351', alpha=0.2, label=layers[3])

ax.legend(loc='lower left')

plt.tight_layout()
plt.savefig('py/Initial continental geotherm.png', dpi=300, transparent=True)
plt.show()
