# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 22:11:09 2020

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

fig, ax1 = plt.subplots(figsize=(28, 20))

    
for a in A:

    Mi = []
    for phi in Phi:
        
        
        Mi.append(1/(np.sqrt(a**2/(1-phi)+1/phi**2)))
    
    
    xy = ( Phi[Mi.index(max(Mi))]+0.03 , max(Mi)+0.004 )
    
    plt.plot(Phi,Mi,color = 'black')
    ax1.annotate("A = %.1f"%(a),xy)


Mi_mf = []


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

##### ##### #####
'''
Alpha2 = []

for phi in Phi:

    Alpha2.append(2*np.arctan(np.sqrt(2*(1-phi)/phi))*180/np.pi)

fig, ax5 = plt.subplots(figsize=(28, 20))

plt.plot(A,Alpha2, color = 'black')

plt.xticks(np.arange(0.0,A_max+0.5,0.5))
plt.xlim(0.0,A_max)

#plt.yticks(np.arange(0.0,1.0+0.02,0.02))
#plt.ylim(0.0,1.0)

plt.ylabel("2α")
plt.xlabel("A")

plt.grid()

plt.savefig('Output/A-Phi.png', dpi= 300, bbox_inches = 'tight')

plt.show()
'''
##### ##### #####

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


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


m_target = 40 #g/s

DeltaP = 5 #bar


Alpha2 = 120 #deg

A_min = 2.0

A_max = 2.0

A_step = 1.0

A_def = [A_min] #np.arange(A_min,A_max+A_step,A_step)


N_in = 4


E = []


for a in A_def:
    
    R_n = -2

    R_n1 = -1
    
    R_in = -1
    
    r_in = -1
    
    Iteration = 1

    It_max = 5
    
    for i in range(len(A)):

        if(abs(a-A[i])<0.01):
    
            Mi_def = Mi_mf[i]

    R_n =  0.475*np.sqrt(m_target*1e-3/(Mi_def*np.sqrt(Rho_p*DeltaP*1e5)))*1000 #mm
                
    
    while(abs(R_n-R_n1)>1e-8 and Iteration <= It_max):
    
        print("Iteração: ",Iteration)
        
        print("E_i = ",abs(R_n-R_n1))
        
        if(R_n1 > 0):
            R_n = R_n1
            a = R_in*R_n/(N_in*r_in**2)
        
        print("\nA = \t\t%.2f"%(a))
        print("\nMi = \t\t%.4f"%(Mi_def))
        
        print("\nR_n = \t\t%.2f \tmm"%(R_n))
        
        print("\nD_n = \t\t%.2f \tmm"%(2*R_n))
        
        C = 1.25
        
        print("\nC = \t\t%.2f"%(C))
        
        R_in = R_n * C
        
        print("\nR_in = \t\t%.2f \tmm"%(R_in))
    
        
        print("\nN_in = \t\t%d"%(N_in))
        
        r_in = np.sqrt(R_in*R_n/(N_in*a))
        
        print("\nr_in = \t\t%.2f \tmm"%(r_in))
        
        print("\nd_in = \t\t%.2f \tmm"%(2*r_in))
        
        l_in = 4*r_in #3~6
        
        print("\nl_in =\t\t %.2f \tmm"%(l_in))
        
        l_n = 2*R_n #0.5~2
        
        print("\nl_n = \t\t%.2f \tmm"%(l_n))
        
        l_s = 3*R_in #>2
        
        print("\nl_s = \t\t%.2f \tmm"%(l_s))
        
        R_s = R_in+r_in
        
        print("\nR_s = \t\t%.2f \tmm"%(R_s))
        
        Re_in = (2/np.pi)*m_target*1e-3/(N_in*r_in*1e-3*Din_visc_p)
        
        print("\nRe_in = \t%d"%(Re_in))
        
        Lambda = 0.3164/(Re_in**0.25)
        
        print("\nLambda = \t%.4f"%(Lambda))
        
        ##### Recalc #####
        
        A_eq = R_in*R_n/(N_in*r_in**2+Lambda/2*R_in*(R_in-R_n))
        
        print("\nA_eq = \t\t%.2f"%(A_eq))
        
        for i in range(len(A)):
        
            if(abs(A_eq-A[i])<0.01):
        
                Mi_eq = Mi_mf[i]
                
                
        print("\nMi_eq = \t%.4f"%(Mi_eq))
        
        for i in range(len(A)):
        
            if(abs(A_eq-A[i])<0.01):
        
                Phi_eq = Phi_mf[i]
        
        Alpha2_eq = np.arctan(np.sqrt(2*(1-Phi_eq)/(Phi_eq)))
        
        print("\nAlpha2_eq = \t%.1f \tdeg"%(2*Alpha2_eq*180/np.pi))
        
        Alpha_in = 90 - 180/np.pi * np.arctan(R_s/l_in)
        
        print("\nAlpha_in = \t%.1f \tdeg"%(Alpha_in))
        
        Ksi_in = -1*Alpha_in/150 + 1.1
        
        print("\nKsi_in = \t%.4f"%(Ksi_in))
        
        Ksi_i = Ksi_in + Lambda*l_in/(2*r_in)
        
        print("\nKsi_i = \t%.4f"%(Ksi_i))
        
        Mi_i = Mi_eq/(np.sqrt(1+Ksi_i*Mi_eq**2*a**2/C**2))
        
        print("\nMi_i = \t\t%.4f"%(Mi_i))
        
        R_n1 = 0.475*np.sqrt(m_target*1e-3/(Mi_i*np.sqrt(Rho_p*DeltaP*1e5)))*1000 #mm
        
        print("\nR_n1 = \t\t%.2f \tmm"%(R_n1))
        
        E.append(abs(R_n-R_n1))
        
        Iteration += 1
        
        
        
        print("\n\n\n")

print("len E: ",len(E))
print("len it: ",len(np.arange(1,Iteration)))


fig, ax4 = plt.subplots(figsize=(14, 7))

plt.plot(np.arange(1,Iteration),E, color = 'black')

#plt.xticks(np.arange(0.0,A_max+0.5,1.0))
#plt.xlim(0.0,A_max)

#plt.yticks(np.arange(0.0,1.0+0.05,0.05))
#plt.ylim(0.0,1.0)

plt.ylabel("E")
plt.xlabel("it.")

plt.grid()

plt.show()

