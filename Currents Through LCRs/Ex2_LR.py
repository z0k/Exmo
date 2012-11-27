from pylab import *
from numpy import *
from scipy.optimize import leastsq

data = loadtxt("Ex2_LR.txt", dtype='float', skiprows=1)

frequency = data[:, 0] * 2. * pi
voltage_inductor = data[:, 1]
voltage_resistor = data[:, 2]
#Convert the phase difference data into radians.
time_delay = data[:, 3]
phase_diff = data[:, 3] * 10 ** -6. * (frequency)
phase_diff_error = phase_diff * sqrt((.20 / time_delay) ** 2 +
                                     (5. * 2. * pi / frequency) ** 2)
angle = ones(len(frequency)) * pi / 2.
#Limit angle between pi and -pi. (No angles are less than -pi,
#so test is not done)
for i in range(0, len(phase_diff)):
    if phase_diff[i] > pi:
        phase_diff[i] -= pi

#Resistances of resistors in ohms.
resistance = [96700., 28740., 473500., 12750., 3114., 517.2]

#Consult report for reasoning behind choice of error values.
voltage_ratio = voltage_inductor / voltage_resistor
#Error propogation for voltage ratio.
ratio = sqrt((0.2 / voltage_inductor) ** 2 +
             (0.2 / voltage_resistor) ** 2) * voltage_inductor / voltage_resistor

#Error propogation for reactance.
inductive_reactance = voltage_ratio * resistance[5]
i_error = sqrt((ratio / (voltage_ratio)) ** 2 +
               (0.05 / resistance[5]) ** 2) * inductive_reactance


#Value of inductance of coil in henries
inductor = 0.04
#Capacitance in Farads
capacitor = 0.022 * 10 ** -6

#Model to compute phase angle based on reactance data.
phi = arctan(inductive_reactance / resistance[5])

#Model to compute reactance.
imp = sqrt(((frequency) * inductor) ** 2 + resistance[5] ** 2)

errorbar((frequency), inductive_reactance, i_error, fmt='r+')
plot((frequency), imp, label='Model')
title(r'LR Circuit Reactance')
plot((frequency), inductive_reactance, linestyle='None', marker='H',
     label='Experimental Data')
legend()
xlabel(r'$\omega$')
ylabel(r'$Z_L$')
show()

errorbar(log(frequency), phase_diff, phase_diff_error, fmt='r+')
title('LR Circuit Phase Angle')
plot(log(frequency), angle, label='Expected Angle for $q_L$')
plot(log(frequency), phi, label=r'Total Phase Angle: $\arctan({\frac{X_L}{R}})$')
plot(log(frequency), phase_diff, linestyle='None', marker='o',
     label='Experimental Data')
legend(loc='lower right')
xlabel(r'$\ln{\omega}$')
ylabel(r'$\phi$')
show()
