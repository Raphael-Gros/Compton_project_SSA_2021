


import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand


me = 9.1E-31
c = 3E8
q = 1.60E-19
e_0 = 8.85E-12
re = 2.817E-15
E0 = 10e6 * q
N = 1000
N_coll = 100
fig,axs1 = plt.subplots(1, 2, subplot_kw = dict(projection = 'polar'))

#fig2, axs2 = plt.subplots(1, 2, constrained_layout=True) 
#plt.figure()

#fig1, axs1 = plt.subplots(2)
#fig1.suptitle('Vertically stacked subplots')

fig2, axs2 = plt.subplots(2)
fig2.suptitle('Vertically stacked subplots')

    
for i in range (N_coll):
    
    if np.any((E0 > 30e3 * q)):  #if Energy smaller then 30e3 stop (asborption)
    
        s = E0 / (me * c**2 + 0.5625 * E0)
        R =rand.rand(N)
        E1 = E0 / (1 + s * R + (2 * E0 / (me * c**2) - s) * R**3)
        theta2 = np.arccos(1 - me * c**2 / E0 * (E0 / E1-1))


        E0 = E1  #reinitialising the energy
        
        #plot-------------------------------------


        axs1[0].hist(theta2,bins = 50, density = True)
        axs1[0].set_title('Theta distribution')
        
        axs1[1].plot(theta2, E1 * 2 * np.pi * np.sin(theta2),'*')
        axs1[1].set_title('Theta distribution (cross-section)',x = 0.7,y = 0.95)
        
        axs2[0].plot(theta2, E1 * 2 * np.pi * np.sin(theta2),'*')
        axs2[0].set_title('Energy vs theta')
        axs2[0].set_xlabel("theta")
        axs2[0].set_ylabel("Energy")
        
        
        axs2[1].hist(theta2, bins = 50, density = True)
        axs2[1].set_title('Theta histogram',x = 0.7,y = 0.95)
        
        #axs2[1,0].plot(E1 * 2 * np.pi * np.sin(theta2), theta2,'*')
        #axs2[1,0].set_title('thata vs Energy')
        #axs2[1,0].set_xlabel("Energy")
        #axs2[1,0].set_ylabel("theta")
