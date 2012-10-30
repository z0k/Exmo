from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import leastsq

data = loadtxt('30degreesweep.txt', dtype='float', skiprows=1)
voltage = data[:, 0]
area = data[:, 1]
current = data[:, 2]
peak_angle = data[:, 3]

#Define constants:
#Thermal coefficient of lightbulb at room temperature.
alpha_0 = 0.0045
#Room temperature in Kelvin.
T_0 = 300.0
#Bulb resistance at room temperature.
R_0 = 1.1
T = zeros(7)
#Initial angle.
Init = 79.7
#Define array for fit parameters.
p = zeros(1)
#Initialize with initial guesses.
p[0] = 5.0
#Array for index of refraction.
n = zeros(7)
#Array for angle data
angle = zeros(7)
#Initial  angle is a given constant. Peak angle was determined by looking at
#the angle at maximum intensity from observations.
angle = Init - peak_angle
angle = angle * pi / 180.
lambd = zeros(7)

#Stefan's constant.
sigma = 5.670400 * 10 ** -8
print sigma

A = 13900.
B = 1.689


def temp(voltage, current):
    return T_0 + ((voltage / current) / R_0 - 1.) / alpha_0


def refrac(angle):
    return sqrt(((2. / sqrt(3.)) * sin(angle) + 1. / 2.) ** 2 + 3. / 4.)


def wavelength(n):
    return sqrt(A / (n - B))


def SBL(T):
    return sigma * T ** 4


#def residuals(sigma, area, T):
#    return area - SBL(T, sigma)


#Generate array for calculated temperature values.
for i in range(0, 7):
    T[i] = temp(voltage[i], current[i])

for i in range(0, 7):
    n[i] = refrac(angle[i])
    lambd[i] = wavelength(n[i])

#p_final, cov_x, info, mesg, success = leastsq(residuals, sigma,
#args=(area, T ** 4), full_output=True)


#print n
#print lambd
wien = mean(lambd * T * 10. ** -9)

print "The value obtained for Wien's displacement constant: ", wien

plot(lambd, n, linestyle='None', marker='o')
show()
print sigma

plot(T, area)
plot(T, SBL(T))

show()
#p_final, cov_x, info, mesg, success = leastsq(residuals, p,
#args=(index_of_refraction, wavelength), full_output=True)

#p_final is a value outputed from the leastsq function.
#y_final = peval(wavelength, p_final)
#Chi-square.
#chi2 = sum((index_of_refraction - y_final) ** 2 / abs(0.3))
#Degrees of freedom.
#dof = len(wavelength) - len(p_final)
#Reduced chi-squared value is calculated.
#print "RMS of residuals (sqrt(chisquared/dof))", sqrt(chi2 / dof)


#for i in range(len(p_final)):
#    print "p[%d] =%12.3f +/- %.4f" % (i, p_final[i], sqrt(cov_x[i, i]))


#title('Something')
#plot(wavelength, index_of_refraction, marker='o', linestyle='None')
#xlabel('Wavelength')
#ylabel('Index of Refraction')
#plot(wavelength, peval(wavelength, p_final))
#show()
