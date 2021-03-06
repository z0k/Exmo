from pylab import *
from numpy import *
from scipy.optimize import leastsq

data = loadtxt("Ex2_LCR.txt", dtype='float', skiprows=1)

frequency = data[:, 0] * 2. * pi
frequency_error = ones(len(frequency)) * 0.5
voltage = data[:, 1]
v_error = ones(len(voltage)) * 0.2
voltage_resistor = data[:, 2]
vr_error = ones(len(voltage_resistor)) * 0.2
time_delay = data[:, 3]
phase_diff = data[:, 3] * 10 ** -6 * frequency


phase_diff_error = phase_diff * sqrt((.20 / time_delay) ** 2 +
                                     (5. * 2. * pi / frequency) ** 2)

#Limit angles between -pi and pi.
for i in range(0, 5):
    phase_diff[i] = phase_diff[i] - pi

#Resistances of resistors in ohms.
resistance = [96700., 28740., 473500., 12750., 3114., 517.2]
resistance_error = 0.05

Z = voltage / voltage_resistor * resistance[5]
ratio = sqrt((v_error / voltage) ** 2 +
             (vr_error / voltage_resistor) ** 2) * voltage / voltage_resistor
Z_error = sqrt((ratio / (voltage / voltage_resistor)) ** 2 +
               (resistance_error / resistance[5]) ** 2) * Z

#Value of inductance of coil in henries
inductor = 0.04
#Capacitance in Farads
capacitor = 0.022 * 10 ** -6

errorbar(log(frequency), Z, Z_error, fmt='r+')
imp = sqrt((frequency * inductor - 1. / (frequency * capacitor)) ** 2 +
            resistance[5] ** 2)
title("LCR circuit")
plot(log(frequency), imp, label='Model')
plot(log(frequency), Z, linestyle='None', marker='o',
     label='Experimental Data')
legend()
xlabel(r'$\ln{\omega}$ (rad s$^{-1}$)')
ylabel(r'$|Z|$ (ohms)')
show()

errorbar(log(frequency), phase_diff, phase_diff_error, fmt='r+')
plot(log(frequency), arctan(((frequency) * inductor - 1. / ((frequency) *
    capacitor)) / resistance[5]), label='Model')
plot(log(frequency), phase_diff, linestyle='None', marker='o',
     label='Experimenta Data')
legend(loc='lower right')
xlabel(r'$\ln{\omega}$')
ylabel(r'$\phi$')
show()
