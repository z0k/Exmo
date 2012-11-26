from pylab import *
from numpy import *
from scipy.optimize import leastsq

data = loadtxt("Ex2_LCR.txt", dtype='float', skiprows=1)

#m_frequency = arange(10,1000000)
frequency = data[:, 0] * 2. * pi
voltage = data[:, 1]
v_error = ones(len(voltage)) * 0.005
voltage_resistor = data[:, 2]
vr_error = ones(len(voltage_resistor)) * 0.005
#Convert the phase angle into degrees.
phase_diff = data[:, 3] * 10 ** -6 * frequency 
#Adjust the phase angle to be within pi.
for i in range(0, 5):
    phase_diff[i] = phase_diff[i] - pi
#Resistances of resistors in ohms.
resistance = [96700., 28740., 473500., 12750., 3114., 517.2]
resistance_error = 0.05
Z = voltage / voltage_resistor * resistance[5]

#Error propogation for uncertainty in impedence measurements.
ratio = sqrt((v_error / voltage) ** 2 + (vr_error / voltage_resistor) ** 2) * voltage / voltage_resistor
Z_error = sqrt((ratio / (voltage / voltage_resistor)) ** 2 + (resistance_error / resistance[5]) ** 2) * Z
  
#Value of inductance of coil in henries
inductor = 0.04
#Capacitance in Farads
capacitor = 0.022 * 10 ** -6
p = zeros(3)
p[0] = 0.
p[1] = 0.
p[2] = 0.


def phi(frequency, p):
    return arctan(((frequency) * inductor - 1. / ((frequency) *
    capacitor)) / (resistance[5]))


def residuals(p, phase_diff, frequency):
    return phase_diff - phi((frequency), p)

        
p_final, cov_x, info, mesg, success = leastsq(residuals, p,
args=(phase_diff, log(frequency)), full_output=True)

#p_final is a value obtained from the leastsq function.
y_final = phi(log(frequency), p_final)
chi2 = sum((Z - y_final) ** 2. / abs(Z_error) ** 2)
#Degrees of freedom.
dof = len(frequency) - len(p_final)
#Reduced chi-squared value is calculated.
print "RMS of residuals (sqrt(chisquared/dof))", sqrt(chi2 / dof)

#for i in range(len(p_final)):
#    print "p[%d] =%12.10e +/- %.11e" % (i, p_final[i], sqrt(cov_x[i, i]))

#errorbar(log(frequency), Z, Z_error, fmt='r+')

title("LCR circuit")

plot(log(frequency), phi(frequency, p_final), label='Model')
plot(log(frequency), phase_diff, linestyle='None', marker='o', label='Experimenta Data')
legend(loc='lower right')
xlabel(r'$\ln{\omega}$')
ylabel(r'$\phi$')
show()
