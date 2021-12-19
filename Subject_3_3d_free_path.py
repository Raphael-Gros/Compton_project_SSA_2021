

### Free path

import numpy as np
import matplotlib.pyplot as plt
import random as rand
from mpl_toolkits.mplot3d import Axes3D

"""Constants"""
N = 2.688E25 #densitÃ© atmospherique
re = 2.81794E-15 #rayon electron
me = 9.1E-31
c=3e8
erest = me*c**2/(1.602E-19) # pour avoir en eV
n = 1 # Number of photons
Nc = 100  # Number of Collisions

"""Functions"""
def CDF(lembda,r): #cumulative distribution function, gives the free path
     return(lembda*np.log(1-r))

def cross_sec(eps,epsp): #cross section per atom
     return(np.pi*re**2*(2*(epsp/eps)**3-4/3*(epsp/eps)**2+2*(epsp/eps)))

def theta(epsp, eps): #angle after collision depending on the energy
     return(np.arccos(1-erest/eps*(eps/epsp-1)))


eps = np.zeros((Nc,n))  # Energy of the photon
epsp = np.zeros((Nc,n))  # Energy after collision
s = np.zeros((Nc,n))

x = np.zeros((Nc,n))
y = np.zeros((Nc,n))
z = np.zeros((Nc,n))

R = np.random.rand(Nc,n)
tita = np.zeros((Nc,n))
fi = np.zeros((Nc,n))
xp = np.random.rand (Nc,n)
yp = np.random.rand (Nc,n)

for j in range (n):
  x[0,j]= 10*np.random.rand()

y[0,:]=0
eps[0,:]=10E6
emin = 30E3


"""3D"""

for i in range(1,Nc):
    for j in range (n):
        if eps[i-1,j]> emin:

            s[i,j] = eps[i-1,j] / (erest + 0.5625*eps[i-1,j])
            eps[i,j] = eps[i-1,j] /(1 + s[i,j]*R[i,j] + (2*eps[i-1,j]/erest-s[i,j])*R[i,j]**3) #energy after collision
            CS = cross_sec(eps[i-1,j],eps[i,j])
            lembda = 1/(CS*N) #mean free path
            cdf = CDF(lembda,R[i,j])
            nu = N*CS*c
            tau = 1/nu
            Pcoll = 1-cdf
            if R[i,j] < Pcoll:
                fi[i,j] = 2*np.pi*R[i,j]
                tita[i,j] = theta(eps[i,j],eps[i-1,j])
                x[i,j] = x[i-1,j] + abs(lembda*np.sin(tita[i,j])*np.cos(fi[i,j]))
                y[i,j] = y[i-1,j] + lembda*np.sin(tita[i,j])*np.sin(fi[i,j])
                z[i,j] = z[i-1,j] + lembda*np.cos(tita[i,j])
        else :
            x[i,j] = np.NaN
            y[i,j] = np.NaN
            z[i,j] = np.NaN

fig = plt.figure()
axes = fig.add_subplot(111, projection='3d')
axes.set_xlabel('X')
axes.set_ylabel('Y')
axes.set_zlabel('Z')
axes.set_zlim(0,10000)
axes.plot_wireframe(x,y,z)

#axes.scatter(x,y,z)
plt.show()

