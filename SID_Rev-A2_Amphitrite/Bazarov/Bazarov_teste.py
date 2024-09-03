# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 22:28:05 2020

@author: A. Goulart
"""

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt

#####  Fluid Properties: Start  #####

Propellant = "LOx"

Rho_p = 1141.0 #kg/m³

Din_visc_p = 2.21e-3 #Pa.s

#####  Fluid Properties: End  #####


A_min = 0.50

A_max1 = 1.50

A_step1 = 0.25

A_max2 = 2.0

A_max3 = 4.0

A_step2 = 1.00

A_max4 = 6.0

A_max5 = 10.00

A_step3 = 4.00


A = np.concatenate([np.arange(A_min,A_max1+A_step1,A_step1), np.arange(A_max2,A_max3+A_step2,A_step2), np.arange(A_max4,A_max5+A_step3,A_step3)])

#print(A)


Phi_min = 0.0
    
Phi_max = 1.0

Phi_step = (Phi_max-Phi_min)/1000

Phi = np.arange(Phi_min,Phi_max,Phi_step)

#fig, ax1 = plt.subplots(figsize=(28, 20))

    
for a in A:

    Mi = []
    for phi in Phi:
        
        
        Mi.append(1/(np.sqrt(a**2/(1-phi)+1/phi**2)))
    
    
    xy = ( Phi[Mi.index(max(Mi))]+0.03 , max(Mi)+0.004 )
    
    #plt.plot(Phi,Mi,color = 'black')
    #ax1.annotate("A = %.1f"%(a),xy)


Mi_mf = []

'''
for phi in Phi:
    
    Mi_mf.append(phi*np.sqrt(phi/(2-phi)))

plt.plot(Phi,Mi_mf, '-.', color = 'black')

plt.xticks(np.arange(0.0,1.0+0.02,0.02))
plt.xlim(0.0,1.0)

plt.yticks(np.arange(0.0,0.65+0.02,0.02))
plt.ylim(0.0,0.65)

plt.ylabel("μ")
plt.xlabel("φ")

plt.grid()

plt.savefig('Output/Phi-Mi.png', dpi= 300, bbox_inches = 'tight')


plt.show()
'''

##### ##### ##### #####

A_min = 0.05

A_max = 25.00

A_step = 0.01

A = np.arange(A_min,A_max+A_step,A_step)

Mi_mf = []

Phi_mf = []

for a in A:

    Mi = []
    for phi in Phi:
        
        
        Mi.append(1/(np.sqrt(a**2/(1-phi)+1/phi**2)))
        
    Mi_mf.append(max(Mi))
    Phi_mf.append(Phi[Mi.index(max(Mi))])
    
'''  
fig, ax2 = plt.subplots(figsize=(28, 20))

plt.plot(A,Mi_mf, color = 'black')

plt.xticks(np.arange(0.0,25.0+2.0,0.5))
plt.xlim(0.0,25.00)

plt.yticks(np.arange(0.0,1.0+0.02,0.02))
plt.ylim(0.0,1.00)

plt.ylabel("μ")
plt.xlabel("A")

plt.grid()

plt.savefig('Output/A-Mi.png', dpi= 300, bbox_inches = 'tight')

plt.show()
'''
##### ##### #####

Alpha2 = []

for phi in Phi_mf:

    Alpha2.append(2*np.arctan(np.sqrt(2*(1-phi)/phi))*180/np.pi)



fig, ax1 = plt.subplots(figsize=(28, 20))

color = '0.0'
ax1.set_xlabel('A')
ax1.set_ylabel('2α [deg]', color=color)
ax1.plot(A, Alpha2, 'k-', markersize=1, label = "2α", color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax1.set_xticks(np.arange(0.0,A_max+0.5,0.5))
ax1.set_xlim(0.0,A_max)

ax1.set_yticks(np.arange(20.0,150.0+10.0,10.0))
ax1.set_ylim(30.0,150.0)

#ax1.grid(color='1.0', linestyle='-', linewidth=3)

#ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = '0.5'
ax2.set_ylabel('μ', color=color)  # we already handled the x-label with ax1
ax2.plot(A, Mi_mf, 'k-.', markersize=1, label = "μ", color=color)
ax2.tick_params(axis='y', labelcolor=color)

ax2.set_xticks(np.arange(0.0,A_max+0.5,0.5))
ax2.set_xlim(0.0,A_max)

ax2.set_yticks(np.arange(0.0,1.0+0.05,0.05))
ax2.set_ylim(0.0,1.00)

#ax2.grid(color='0.5', linestyle='--', linewidth=1)

'''
l = ax1.get_ylim()
l2 = ax2.get_ylim()
f = lambda x : l2[0]+(x-l[0])/(l[1]-l[0])*(l2[1]-l2[0])
ticks = f(ax1.get_yticks())
ax2.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
'''


for ax in [ax1, ax2]:
    ax.xaxis.grid(True, which='both')
    ax.yaxis.grid(True, which='both')

    
fig.legend(loc = 'upper center')
fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.savefig('Output/A-Alpha-Mi.png', dpi= 300, bbox_inches = 'tight')

plt.show()

'''
fig, ax1 = plt.subplots(figsize=(28, 20))

ax1.plot(A,Alpha2, color = 'black')

#ax1.xticks(np.arange(0.0,A_max+0.5,0.5))
#ax1.xlim(0.0,A_max)

#ax1.yticks(np.arange(20.0,150.0+5.0,5.0))
#ax1.ylim(30.0,150.0)


ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylabel("2α [deg]")
ax1.set_xlabel("A")

ax2 = ax1.twinx()

ax2.plot(A,Mi_mf, color = 'black')
ax2.set_ylabel("μ")
ax2.tick_params(axis='y', labelcolor=color)


#plt.grid()

#plt.savefig('Output/A-Phi.png', dpi= 300, bbox_inches = 'tight')

fig.tight_layout()

plt.show()
'''
##### ##### #####
'''
fig, ax3 = plt.subplots(figsize=(28, 20))

plt.plot(A,Phi_mf, color = 'black')

plt.xticks(np.arange(0.0,A_max+0.5,0.5))
plt.xlim(0.0,A_max)

plt.yticks(np.arange(0.0,1.0+0.02,0.02))
plt.ylim(0.0,1.0)

plt.ylabel("φ")
plt.xlabel("A")

plt.grid()

plt.savefig('Output/A-Phi.png', dpi= 300, bbox_inches = 'tight')

plt.show()
'''