### Free path
import numpy as np
import matplotlib.pyplot as plt
import random as rand

"""Constants"""
N = 2.688E25 #densité atmospherique
re = 2.81794E-15 #rayon electron
me = 9.1E-31
c=3e8
h = 6.626E-34
erest = me*c**2/(1.602E-19) # pour avoir en eV
def CDF(lembda,r): #cumulative distribution function, gives the free path
     return(lembda*np.log(1-r))

def cross_sec(eps,epsp): #cross section per atom
     return(np.pi*re**2*(2*(epsp/eps)**3-4/3*(epsp/eps)**2+2*(epsp/eps)))

def theta(epsp, eps): #angle after collision depending on the energy
     return(np.arccos(1-erest/eps*(eps/epsp-1)))



n = 200000 # Number of photons
Nc = 30  # Number of steps

eps = np.zeros((Nc,n))  # Energy of the photon  
epsp = np.zeros((Nc,n))  # Energy after collision
s = np.zeros((Nc,n)) 
x = np.zeros((Nc,n)) 
y = np.zeros((Nc,n)) 
R = np.random.rand (Nc,n)
tita = np.zeros((Nc,n)) 
t = np.zeros((Nc,n))
#yp = np.random.rand (Nc,n)
counter = 0
espec = np.zeros(n)

for j in range (n):
  x[0,j]= 10*np.random.rand()+10000
  
y[0,:]=0
eps[0,:]=10E6
v = 30E3

for j in range (n):
  for i in range(1,Nc):

      
          if eps[i-1,j]> v:

              s[i,j] = eps[i-1,j] / (erest + 0.5625*eps[i-1,j])
              eps[i,j] = eps[i-1,j] /(1 + s[i,j]*R[i,j] + (2*eps[i-1,j]/erest-s[i,j])*R[i,j]**3) #energy after collision
              CS = cross_sec(eps[i-1,j],eps[i,j])
              lembda = 1/(CS*N) #mean free path
              cdf = CDF(lembda,R[i,j])
              nu = N*CS*c
              tau = 1/nu
              Pcoll = 1-cdf
              t[i,j] = t[i-1,j] + lembda/c
              if R[i,j] < Pcoll:
                  tita[i,j] = theta(eps[i,j],eps[i-1,j])
                  #print(lembda*np.cos(tita[i,j]))
                  x[i,j] = x[i-1,j] + lembda*np.cos(tita[i,j])
                  y[i,j] = y[i-1,j] + lembda*np.sin(tita[i,j])
                  
              else :
                  x[i,j] = np.NaN
                  y[i,j] = np.NaN
          if  y[i,j] >= 15000 and y[i,j]<= 16000 :
              if x[i,j] <= 10000 and x[i,j] >= 0 :
                 espec[j] = eps [i,j] 
                 counter = counter + 1      
espec = espec[espec!=0]             
freq = espec*1.602E-19 / h  
#plt.plot(x,y,'-')               
#plt.plot(freq,espec,'.')
print('Nombre de photon par seconde',counter)            
histvals=plt.hist(espec,bins=100)
histx,histy=histvals[0],histvals[1]
#plt.plot(histy,histx,'.')
plt.show()

