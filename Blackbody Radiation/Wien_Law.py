from numpy import *
from scipy import *
from pylab import *
from scipy.optimize import leastsq

data = loadtxt('30degreesweep.txt', dtype='float', skiprows=1)
voltage = data[:, 0]
voltage_error = ones(7) * 0.0005
area = data[:, 1]
current = data[:, 2]
current_error = ones(7) * 0.0005
peak_angle = data[:, 3]

#Define constants:
#Thermal coefficient of lightbulb at room temperature.
alpha_0 = 0.0045
#Room temperature in Kelvin.
T_0 = 300.0
#Bulb resistance at room temperature.
R_0 = 1.1
T = zeros(7)
ratio = sqrt( (voltage_error / voltage) ** 2 + (current_error / current) ** 2) * voltage / current
#Initial angle.
Init = 79.7
#Define array for fit parameter for Wien's Law.
p = zeros(1)
#Initialize with initial guesses.
p[0] = 0.002
#Array for index of refraction.
n = zeros(7)
#Array for angle data
angle = zeros(7)
#Initial  angle is a given constant. Peak angle was determined by looking at
#the angle at maximum intensity from observations.
angle = Init - peak_angle
angle = angle * pi / 180.
lambd = zeros(7)
lambd_error = zeros(7)
angle_error = (Init - array([0.25870, 0.14490, 0.16557, 0.31044, 0.23800, 0.26906, 0.25870])) * pi / 180.


#Stefan's constant.
sigma = 5.670400 * 10 ** -8

A = 13900.
B = 1.689


def temp(voltage, current):
    return T_0 + ((voltage / current) / R_0 - 1.) / alpha_0


def temp_error(ratio):
    return T_0 + ((ratio) / R_0 - 1.) / alpha_0


def refrac(angle):
    return sqrt(((2. / sqrt(3.)) * sin(angle) + 1. / 2.) ** 2 + 3. / 4.)


def wavelength(n):
    return sqrt(A / (n - B))


def SBL(T):
    return sigma * T ** 4


def wien_law(T, p):
    return p[0] * 1. / T


def residuals(p, lambd, T):
    return lambd - wien_law(T, p)


T_error = temp_error(ratio)

#Generate array for calculated temperature values.
for i in range(0, 7):
    T[i] = temp(voltage[i], current[i])

for i in range(0, 7):
    n[i] = refrac(angle[i])
    lambd[i] = wavelength(n[i])
    #The angle error must be converted to wavelength error in nanometers. To
    #do this, it is passed to the refrac function, which in turn is passed
    #to the wavelength function.
    lambd_error[i] = wavelength(refrac(angle_error[i]))


p_final, cov_x, info, mesg, success = leastsq(residuals, p,
args=(lambd, T), full_output=True)

errorbar(T, lambd, lambd_error, fmt='r+')

xlabel('Temperature (Kelvin)')
ylabel('Wavelength (nanometers)')
plot(T, lambd, linestyle="None", marker='o')
plot(T, wien_law(T, p_final))
show()

wien = mean(lambd * T * 10. ** -9)

print "The value obtained for Wien's displacement constant by taking the average: %e" % wien, "\n"


#p_final is a value outputed from the leastsq function.
y_final = wien_law(lambd, p_final)
#Chi-square.
chi2 = sum((lambd - y_final) ** 2 / (abs(lambd_error) ** 2))
#Degrees of freedom.
dof = len(lambd) - len(p_final)
#Reduced chi-squared value is calculated.
print "RMS of residuals (sqrt(chisquared/dof))", sqrt(chi2 / dof)


for i in range(len(p_final)):
    print "p[%d] =%12.3f +/- %.4f" % (i, p_final[i], sqrt(cov_x[i, i]))

print "\nThe value obtained for Wien's displacement constant given by p_final: %e" % (p_final[0] * 10. ** -9)
