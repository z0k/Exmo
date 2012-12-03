from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import leastsq

data = loadtxt('data.txt', dtype='float', skiprows=1)
#Frequency in Hz.
frequency = data[:, 0] * 10 ** 6 
angle = data[:, 1] * pi / 180.

#Wavelength of the sodium lamp.
dL = 1. * 2.54 * 10 ** -2 / 2500
lambda_L = zeros(2)
lambda_L[0] = dL * sin(3.4 * pi / 180.) / 1.
lambda_L[1] = dL * sin(6.9 * pi / 180.) / 2.

#Order of diffraction.
m = 1.

lambda_L_avg = (lambda_L[0] + lambda_L[1]) / 2.
#Angle error in radians.
angle_error = 0.016 * pi / 180.
lambda_L_error = dL * sin(angle_error)
#This is the wavelength of sound under water.
spacing = zeros (len(angle))
spacing = m * lambda_L_avg / sin(angle)
spacing_error = zeros(len(angle))
spacing_error = m * sqrt( (lambda_L_error / lambda_L_avg) ** 2 + (angle_error / angle) ** 2 ) * spacing

p = zeros(1)
p[0] = 300.


def peval(frequency, p):
    return p[0] / frequency


def residuals(p, spacing, freqeuncy):
    return spacing - peval(frequency, p)


#Least squares.
p_final, cov_x, info, mesg, success = leastsq(residuals, p,
args=(spacing, frequency), full_output=True)

 
#p_final is a value obtained from the leastsq function.
y_final = peval(frequency, p_final)
#Chi-square.
chi2 = sum((spacing - y_final) ** 2 / (abs(spacing_error) ** 2))
#Degrees of freedom.
dof = len(spacing) - len(p_final)
#Reduced chi-squared value is calculated.
print "RMS of residuals (sqrt(chisquared/dof))", sqrt(chi2 / dof)


for i in range(len(p_final)):
    print "p[%d] =%12.3f +/- %.4f" % (i, p_final[i], sqrt(cov_x[i, i]))


v_sound = p_final[0]
v_sound_error = sqrt(cov_x[0, 0])

#Density of water.
rho = 1000.
#Bulk modulus of water.
b_modulus = v_sound ** 2 * rho
#Bulk modulus error.
b_modulus_error = (2. * v_sound_error * v_sound_error) * rho

print "\nThe velocity of sound in water is %.3f +/- %.4f meters per second." % (v_sound, v_sound_error)
print "The bulk modulus of water is %.5e +/- %.6e pascals." % (b_modulus, b_modulus_error)

errorbar(frequency, spacing, spacing_error, fmt='r+')
plot(frequency, peval(frequency, p_final))
plot(frequency, spacing, linestyle='None', marker='o')
xlabel("frequency (Hz)")
ylabel("wavelength of sound (m)")
show()
