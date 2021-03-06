from pylab import *
from numpy import *
from scipy.optimize import leastsq

data = loadtxt("Ex2_LR.txt", dtype='float', skiprows=1)


#Frequency in Hertz converted into angular frequency.
frequency = data[:, 0] * 2. * pi
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
reactance_inductor = frequency * inductor
reactance_capacitor = 1. / (frequency * capacitor)

def impedance(reactance, resistance):
    return (reactance ** 2 + resistance ** 2) ** 0.5


plot(log(frequency), log(impedance(reactance_inductor, resistance[5])))
plot(log(frequency), log(voltage_ratio * resistance[5]))
show()
