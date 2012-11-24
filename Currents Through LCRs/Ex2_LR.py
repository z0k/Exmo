from pylab import *
from numpy import *
from scipy.optimize import leastsq

data = loadtxt("Ex2_LR.txt", dtype='float', skiprows=1)

frequency = data[:, 0]
voltage_inductor = data[:, 1]
voltage_resistor = data[:, 2]
#Convert the phase difference data into seconds.
phase_diff = data[:, 3] * 10 ** -3.

voltage_ratio = voltage_inductor / voltage_resistor

#Resistances of resistors in ohms.
resistance = [96700., 28740., 473500., 12750., 3114., 517.2]
#Value of inductance of coil in henries
inductor = 0.3
#Capacitance in Farads
capacitor = 0.022 * 10 ** -6
reactance_inductor = 2 * pi * frequency * inductor
reactance_capacitor = 1. / (2 * pi * frequency * capacitor)


def impedance(reactance, resistance):
    return sqrt(reactance ** 2 + resistance ** 2)


#plot(frequency, impedance(reactance_inductor, resistance[5]))
title(r'LR Circuit Experimental data')
plot(frequency, voltage_ratio * resistance[5], linestyle='None', marker='o')
xlabel(r'$f$')
ylabel(r'$Z$')
show()
