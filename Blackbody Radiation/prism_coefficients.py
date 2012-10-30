from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import leastsq

data = loadtxt('Refraction_wavelength.txt', dtype='float', skiprows=1)
index_of_refraction = data[:, 0]
wavelength = data[:, 1]

#Define array for fit parameters.
p = zeros(2)
#Initialize with initial guesses.
p[0] = 500.0
p[1] = 2.0


def peval(wavelength, p):
    return p[0] / wavelength ** 2 + p[1]


def residuals(p, index_of_refraction, wavelength):
    return index_of_refraction - peval(wavelength, p)


p_final, cov_x, info, mesg, success = leastsq(residuals, p,
args=(index_of_refraction, wavelength), full_output=True)

#p_final is a value outputed from the leastsq function.
y_final = peval(wavelength, p_final)
#Chi-square.
chi2 = sum((index_of_refraction - y_final) ** 2 / abs(0.3))
#Degrees of freedom.
dof = len(wavelength) - len(p_final)
#Reduced chi-squared value is calculated.
print "RMS of residuals (sqrt(chisquared/dof))", sqrt(chi2 / dof)


for i in range(len(p_final)):
    print "p[%d] =%12.3f +/- %.4f" % (i, p_final[i], sqrt(cov_x[i, i]))


title('Something')
plot(wavelength, index_of_refraction, marker='o', linestyle='None')
xlabel('Wavelength')
ylabel('Index of Refraction')
plot(wavelength, peval(wavelength, p_final))
show()
