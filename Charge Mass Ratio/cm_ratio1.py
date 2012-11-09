from numpy import *
from pylab import *
from scipy.optimize import leastsq

data = loadtxt("charge_mass_ratio.txt", dtype='float', skiprows=1)

voltage = data[:5, 0]
current = data[:5, 1]
diameter = data[:5, 2]

voltage_error = ones(len(voltage)) * 0.05
current_error = ones(len(current)) * 0.0005
diameter_error = ones(len(diameter)) * 0.5
voltage_root_error = 0.5 * (voltage) ** (-.5) * 0.05

#Convert to radius and convert to meters.
radius = diameter / 200.
radius_error = diameter_error / 200.
#Y value is computed as this expression to check linear relationship with
#a changing current as the independent variable.
y = sqrt(voltage) / radius
#Error propogation in expression for y.
y_error = (sqrt(voltage) / radius) * sqrt((voltage_root_error / voltage) ** 2 + (0.0025 / radius) ** 2)

#Fit parameters initialized with guess of 1.
p = zeros(2)
p[0] = 1.
p[1] = 1.

#Constants; mu is described in the document. n is the number of turns in
#each coil, and R is the radius of the coils.
mu = 4. * pi * 10. ** -7
n = 130.
R = 0.161
#Expression for the constant k in equation 8 from document.
k = (mu * n / R) * (1 / sqrt(2.)) * (4. / 5.) ** (3. / 2.)


#The function representing equation 8 from the document.
def eight(current, p):
    return p[0] * k * (current - p[1])


def residuals(p, y, current):
    return y - eight(current, p)


p_final, cov_x, info, mesg, success = leastsq(residuals, p,
args=(y, current), full_output=True)

y_final = eight(current, p_final)
chi2 = sum((y - y_final) ** 2 / abs(y_error ** 2))
#Degrees of freedom.
dof = len(current) - len(p_final)
#Reduced chi-squared value is calculated.
print "RMS of residuals (sqrt(chisquared/dof))", sqrt(chi2 / dof)
#Output of fit statistics.
for i in range(len(p_final)):
    print "p[%d] =%12.4e +/- %.5e" % (i, p_final[i], sqrt(cov_x[i, i]))

errorbar(current, y, y_error, fmt='r+')

#p[0] is squared in order to obtain the charge mass ratio.
em_ratio = p_final[0] ** 2
#Error propogation for charge mass ratio.
em_ratio_error = 2. * p[0] * sqrt(cov_x[0, 0])
print "The calculated value of (e/m) is%12.4e +/_ %.5e." % (em_ratio, em_ratio_error)

suptitle("Data Set 1")
title("Charge Mass Ratio of an Electron")
plot(current, y, linestyle='None', marker='o')
plot(current, eight(current, p_final))
ylabel(r'$\frac{\sqrt{V}}{r}$', fontsize=18)
xlabel(r'$I$', fontsize=14)
show()
