from pylab import *
from numpy import *
x=linspace(-pi,pi,1000)
n=20
y=zeros(1000)
for i in x:
    temp=0
    for j in range(1,20):
        temp+=cos(2*j*i)/(2*j**2-1)
    y[i]=(2/pi)-(4/pi)*temp
#for i in x:
#    temp=0
#    for j in range(1,20):
#        temp+=-(4/pi)*(cos(j*i)*cos(j*pi)/(j^2-1))
#    y[i]=(2/pi)+temp
plot(x,y)
show()

plot(x, cos(x))
show()

