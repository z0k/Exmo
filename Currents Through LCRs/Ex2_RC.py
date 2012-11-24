from pylab import *
from numpy import *
from scipy.optimize import leastsq

data = loadtxt("Ex2_RC.txt", dtype='float', skiprows=1)

<<<<<<< HEAD
#Frequency in Hertz converted into angular frequency.
frequency = data[:, 0]  * 2 * pi
=======
frequency = data[:, 0] * 2. * pi
>>>>>>> 0e852ebf64571be7aeb2786ec75882776c3b9c1e
voltage_capacitor = data[:, 1]
vc_error = ones(len(voltage_capacitor)) * 0.005
voltage_resistor = data[:, 2]
<<<<<<< HEAD

voltage_ratio = voltage_capacitor / voltage_resistor

=======
vr_error = ones(len(voltage_resistor)) * 0.005
#Convert the phase difference data into seconds.
phase_diff = data[:, 3] * 10 ** -6.
>>>>>>> 0e852ebf64571be7aeb2786ec75882776c3b9c1e
#Resistances of resistors in ohms.
resistance = [96700., 28740., 473500., 12750., 3114., 517.2]
resistance_error = 0.05
Z = resistance[5] * voltage_capacitor / voltage_resistor

ratio = sqrt((vc_error / voltage_capacitor) ** 2 + (vr_error / voltage_resistor) ** 2) * voltage_capacitor / voltage_resistor
Z_error = sqrt((ratio / (voltage_capacitor / voltage_resistor)) ** 2 + (resistance_error / resistance[5]) ** 2) * Z

p = zeros(2)
p[0] = 0.5
p[1] = 1.

#Value of inductance of coil in henries
inductor = 0.45
#Capacitance in Farads
capacitor = 0.022 * 10 ** -6
reactance_inductor = log(frequency) * inductor
reactance_capacitor = 1. / (log(frequency) * capacitor)

<<<<<<< HEAD
phase_diff = data[3:-2, 3] / frequency[3:-2]
print phase_diff


def impedance(reactance, resistance):
    return (reactance ** 2 + resistance ** 2) ** 0.5
=======
>>>>>>> 0e852ebf64571be7aeb2786ec75882776c3b9c1e

def impedance(reactance, p):
    return p[1] * (reactance ** 2 + resistance[5] ** 2) ** p[0]


<<<<<<< HEAD
plot((frequency[3:-2]), (phase_diff))
show()

=======
def residuals(p, Z, reactance):
    return Z - impedance(reactance, p)


p_final, cov_x, info, mesg, success = leastsq(residuals, p,
args=(Z, reactance_capacitor), full_output=True)

#p_final is a value outputed from the leastsq function.
y_final = impedance(reactance_capacitor, p_final)
#Chi-square. IGNORE!! Uncertainty was not given in question.
chi2 = sum((Z - y_final) ** 2 / abs(0.3))
#Degrees of freedom.
dof = len(frequency) - len(p_final)
#Reduced chi-squared value is calculated.
print "RMS of residuals (sqrt(chisquared/dof))", sqrt(chi2 / dof)

#for i in range(len(p_final)):
#    print "p[%d] =%12.10f +/- %.11f" % (i, p_final[i], sqrt(cov_x[i, i]))

errorbar(log(frequency), Z, Z_error, fmt='r+')

title("Semi-log coordinates of RC circuit")
plot(log(frequency), impedance(reactance_capacitor, p_final), linestyle='-', marker='x')
plot(log(frequency), Z, linestyle='None', marker='.')
xlabel(r'$f$')
ylabel(r'$Z$')
show()
>>>>>>> 0e852ebf64571be7aeb2786ec75882776c3b9c1e
