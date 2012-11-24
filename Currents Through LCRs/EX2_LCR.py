from pylab import *
from numpy import *
from scipy.optimize import leastsq

data = loadtxt("Ex2_LCR.txt", dtype='float', skiprows=1)

frequency = data[:, 0]
voltage = data[:, 1]
voltage_resistor = data[:, 2]
#Convert the phase difference data into seconds.
phase_diff = data[:, 3] * 10 ** -6.

voltage_ratio = voltage / voltage_resistor

#Resistances of resistors in ohms.
resistance = [96700., 28740., 473500., 12750., 3114., 517.2]
#Value of inductance of coil in henries
inductor = 0.3
#Capacitance in Farads
capacitor = 0.022 * 10 ** -6
reactance_inductor = 2 * pi * frequency * inductor
reactance_capacitor = 1. / (2 * pi * (frequency) * capacitor)
reactance = reactance_inductor - reactance_capacitor


def impedance(reactance, resistance):
    return sqrt(reactance ** 2 + resistance ** 2)


title("None")
#plot(frequency, impedance(reactance, resistance[5]), linestyle='None', marker='o')
plot(frequency, voltage_ratio * resistance[5], linestyle='None', marker='H', markersize=10)
xlabel(r'$f$')
ylabel(r'$Z$')
show()
