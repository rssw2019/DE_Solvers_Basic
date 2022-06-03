#This program solves the generalized KdV Equations
#Codes written by Ramjee Sharma
#All rights reserved
#6/3/2022
#Import Python libraries
import numpy as np
import matplotlib.pyplot as plt
import time

#Define global parameters
a=0
b=1
c=0
d=0
nSteps=500
dt=0.001

#space discritization
N=256
x0=0
xN=2*np.pi 
x=np.linspace(x0,xN,N, endpoint=False)
u=np.sin(x)

#Wave Numbers
k=np.fft.fftfreq(int(N))*int(N)
#1d RK4 Routine,
def rk4_1d(x,dt):
    k1=dt*f(x)
    k2=dt*f(x+0.5*k1)
    k3=dt*f(x+0.5*k2)
    k4=dt*f(x+k3)
    x=x+(1/6)*(k1+2*k2+2*k3+k4)#Time integration
    return x

#The model u_t=F(a,b,c,d,u)
#u,ux,uxx,uxxx in physycal space
#u1,u1x,u1xx,u1xxx in Fourier Space space
def f(u):
    u1=np.fft.fft(u) #to Fourier Space
    u1x=1j*k*u1
    u1xx=1j*k*u1x 
    u1xxx=1j*k*u1xx 
    ux=np.real(np.fft.ifft(u1x))
    uxx=np.real(np.fft.ifft(u1xx))
    uxxx=np.real(np.fft.ifft(u1xxx))
    return -(a*ux+2*b*u*ux+c*uxxx-d*uxx)

#Time integration using 4th order RK method
for i in range (nSteps):
    u=rk4_1d(u,dt)

#Plotting the results
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(x, u, label='x-u')
plt.title('evolution of pde')
ax.legend()
plt.show()
