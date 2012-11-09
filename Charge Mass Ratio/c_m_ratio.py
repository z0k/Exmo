from numpy import *
from pylab import *
from scipy.optimize import leastsq

data = loadtxt("charge_mass_ratio.txt", dtype='float', skiprows=1)

voltage1 = data[:5, 0]
current1 = data[:5, 1]
diameter1 = data[:5, 2]

voltage2 = data[6:12, 0]
current2 = data[6:12, 1]
diameter2 = data[6:12, 2]

voltage3 = data[13:, 0]
current3 = data[13:, 1]
diameter3 = data[13:, 2]

radius1 = diameter1 / 2.
radius2 = diameter2 / 2.
radius3 = diameter3 / 3.

y1 = sqrt(voltage1) / radius1
y2 = sqrt(voltage2) / radius2
y3 = sqrt(voltage3) / radius3

p1 = zeros(2)
p1[0] = 1.
p1[1] = 1.
p2 = zeros(2)
p2[0] = 1.
p2[1] = 1.
p3 = zeros(2)
p3[0] = 1.
p3[1] = 1.


#Constants.
mu = 4. * pi * 10. ** -7
n = 130.
R = 32.2

k = 1. / (mu * n / R) * sqrt(2.) * (4. / 5.) ** (3. / 2.)


def eight(current, p):
    return p[0] * k * (current - p[1])


#def eight(current, p):
#    return current


def residuals(p, y, current):
    return y - eight(current, p)


p_final1, cov_x1, info, mesg, success = leastsq(residuals, p1,
args=(y1, current1), full_output=True)
p_final2, cov_x2, info, mesg, success = leastsq(residuals, p2,
args=(y2, current2), full_output=True)
p_final3, cov_x3, info, mesg, success = leastsq(residuals, p3,
args=(y3, current3), full_output=True)

plot(y1, current1)
plot(y1, eight(current1, p_final1))
show()

plot(y2, current2)
plot(y2, eight(current2, p_final2))
show()

plot(y3, current3)
plot(y3, eight(current3, p_final3))
show()
