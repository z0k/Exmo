from pylab import *
from numpy import *
from scipy.optimize import leastsq

data = loadtxt("Ex2_RC.txt", dtype='float', skiprows=1)

frequency = data[3:-2, 0] * 2. * pi
voltage_capacitor = data[3:-2, 1]
vc_error = 0.2
voltage_resistor = data[3:-2, 2]
vr_error = 0.2

#Convert the phase difference data into radians.
time_delay = data[3:-2, 3]
phase_diff = data[3:-2, 3] * 10 ** -6. * (frequency) - pi
phase_diff_error = phase_diff * sqrt((.20 / time_delay) ** 2 +
                                     (5. * 2. * pi / frequency) ** 2)
angle = ones(len(frequency)) * -pi / 2.

#Resistances of resistors in ohms.
resistance = [96700., 28740., 473500., 12750., 3114., 517.2]
resistance_error = 0.05
Z = resistance[5] * voltage_capacitor / voltage_resistor

ratio = sqrt((vc_error / voltage_capacitor) ** 2 +
             (vr_error / voltage_resistor) ** 2) * voltage_capacitor / voltage_resistor
Z_error = sqrt((ratio / (voltage_capacitor / voltage_resistor)) ** 2
               + (resistance_error / resistance[5]) ** 2) * Z

errorbar(log(frequency), Z, Z_error, fmt='r+')

imp = sqrt(1. / ((frequency) * 0.022 * 10 ** -6) ** 2 + resistance[5] ** 2)

title("RC circuit Reactance")
plot(log(frequency), imp, linestyle='-', label='Model')
plot(log(frequency), Z, linestyle='None', marker='.',
     label='Experimental Data')
legend()
xlabel(r'$\ln{\omega}$')
ylabel(r'$|Z|$')
show()

#Model to compute phase angle based on reactance data.
phi = arctan(Z / resistance[5])

errorbar(log(frequency), phase_diff, phase_diff_error, fmt='r+')
title('RC Circuit Phase Angle')
plot(log(frequency), angle, label='Expected Angle for $q_R$')
plot(log(frequency), phi,
     label=r'Total phase angle: $\arctan({\frac{X_R}{R}})$')
plot(log(frequency), phase_diff, linestyle='None', marker='o',
     label='Experimental Data')
legend(loc='upper right')
xlabel(r'$\ln{\omega}$')
ylabel(r'$\phi$')
show()
