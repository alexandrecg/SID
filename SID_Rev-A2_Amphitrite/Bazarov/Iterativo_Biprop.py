# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 22:11:09 2020

@author: A. Goulart
"""

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt



#####     #####     #####     #####     #####     #####     #####     #####
#####     #####     #####     Estágio de Oxidante           #####     #####
#####     #####     #####     #####     #####     #####     #####     #####



#####  Fluid Properties: Start  #####

Propellant = "LOx"



if(Propellant == "LOx"):

    Rho_p = 1141.0 #kg/m³
    
    Din_visc_p = 2.21e-3 #Pa.s
    
    m_target = 40 #g/s


elif(Propellant == "Ethanol 95%"):


    Rho_p = 800.0 #kg/m³
    
    Din_visc_p = 0.001095 #Pa.s

    m_target = 25 #g/s

elif(Propellant == "H2O 20ºC"):
    
    Rho_p = 965.31  #kg/m³

    Din_visc_p = 0.0008921 #Pa.s
    
    m_target = 50 #g/s



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

#fig, ax1 = plt.subplots(figsize=(14, 7))

    
for a in A:

    Mi = []
    for phi in Phi:
        
        
        Mi.append(1/(np.sqrt(a**2/(1-phi)+1/phi**2)))
    
    
    xy = ( Phi[Mi.index(max(Mi))]+0.03 , max(Mi)+0.004 )
    
    #plt.plot(Phi,Mi,color = 'black')
    #ax1.annotate("A = %.1f"%(a),xy)


Mi_mf = []


for phi in Phi:
    
    Mi_mf.append(phi*np.sqrt(phi/(2-phi)))

'''
plt.plot(Phi,Mi_mf, '-.', color = 'black')

plt.xticks(np.arange(0.0,1.0+0.1,0.05))
plt.xlim(0.0,1.0)

plt.yticks(np.arange(0.0,0.5+0.15,0.05))
plt.ylim(0.0,0.65)

plt.ylabel("μ")
plt.xlabel("φ")

plt.grid()

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
fig, ax2 = plt.subplots(figsize=(14, 7))

plt.plot(A,Mi_mf, color = 'black')

plt.xticks(np.arange(0.0,25.0+2.0,1.0))
plt.xlim(0.0,25.00)

plt.yticks(np.arange(0.0,0.8+0.15,0.05))
plt.ylim(0.0,0.65)

plt.ylabel("μ")
plt.xlabel("A")

plt.grid()

plt.show()
'''


'''
fig, ax3 = plt.subplots(figsize=(14, 7))

plt.plot(A,Phi_mf, color = 'black')

plt.xticks(np.arange(0.0,A_max+0.5,1.0))
plt.xlim(0.0,A_max)

plt.yticks(np.arange(0.0,1.0+0.05,0.05))
plt.ylim(0.0,1.0)

plt.ylabel("φ")
plt.xlabel("A")

plt.grid()

plt.show()
'''


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


DeltaP = 5 #bar


Alpha2 = 120 #deg

A_min = 6.0

A_max = 6.0

A_step = 1.0

A_def = [A_min] #np.arange(A_min,A_max+A_step,A_step)


C = 1.50



N_in = 4


E = []


for a in A_def:
    
    R_n = -2

    R_n1 = -1
    
    R_in = -1
    
    r_in = -1
    
    Iteration = 1

    It_max = 1000
    
    for i in range(len(A)):

        if(abs(a-A[i])<0.01):
    
            Mi_def = Mi_mf[i]

    R_n =  0.475*np.sqrt(m_target*1e-3/(Mi_def*np.sqrt(Rho_p*DeltaP*1e5)))*1000 #mm
                
    
    while(abs(R_n-R_n1)>1e-8 and Iteration <= It_max):
    
        
        
        if(R_n1 > 0):
            R_n = R_n1
            a = R_in*R_n/(N_in*r_in**2)
        
        
        R_in = R_n * C

        
        r_in = np.sqrt(R_in*R_n/(N_in*a))
        
        
        
        l_in = 4*r_in #3~6
        
        
        
        l_n = 2*R_n #0.5~2
        
        
        
        l_s = 3*R_in #>2
        
        
        
        R_s = R_in+r_in
        
        
        
        Re_in = (2/np.pi)*m_target*1e-3/(N_in*r_in*1e-3*Din_visc_p)
        
        
        
        Lambda = 0.3164/(Re_in**0.25)
        
        
        
        ##### Recalc #####
        
        A_eq = R_in*R_n/(N_in*r_in**2+Lambda/2*R_in*(R_in-R_n))
        
        
        
        for i in range(len(A)):
        
            if(abs(A_eq-A[i])<0.01):
        
                Mi_eq = Mi_mf[i]
                
                
        
        
        for i in range(len(A)):
        
            if(abs(A_eq-A[i])<0.01):
        
                Phi_eq = Phi_mf[i]
        
        Alpha2_eq = np.arctan(np.sqrt(2*(1-Phi_eq)/(Phi_eq)))
        
        
        
        Alpha_in = 90 - 180/np.pi * np.arctan(R_s/l_in)
        
        
        
        Ksi_in = -1*Alpha_in/150 + 1.1
        
        
        
        Ksi_i = Ksi_in + Lambda*l_in/(2*r_in)
        
        
        
        Mi_i = Mi_eq/(np.sqrt(1+Ksi_i*Mi_eq**2*a**2/C**2))
        
        
        
        R_n1 = 0.475*np.sqrt(m_target*1e-3/(Mi_i*np.sqrt(Rho_p*DeltaP*1e5)))*1000 #mm
        
        I = 23.7/((R_in/R_n)**2.7*N_in**1.34*Phi_eq**1.1*(l_s/(2*R_s))**0.15)
        
        
        E.append(abs(R_n-R_n1))
        
        Iteration += 1
        
        
     
        
        



print('\033[4m' + '\t Estágio de Oxidante \t\n\n' + '\033[0m')
        

print("Iteração: ",Iteration-1)

print("E_i = ",abs(R_n-R_n1))



print("\nInjetor Ideal:")

print("\nA = \t\t%.2f"%(a))

print("\nMi = \t\t%.4f"%(Mi_def))

print("\nR_n = \t\t%.2f \tmm"%(R_n))

print("\nD_n = \t\t%.2f \tmm"%(2*R_n))

print("\nC = \t\t%.2f"%(C))

print("\nR_in = \t\t%.2f \tmm"%(R_in))

print("\nN_in = \t\t%d"%(N_in))

print("\nr_in = \t\t%.2f \tmm"%(r_in))

print("\nd_in = \t\t%.2f \tmm"%(2*r_in))

print("\nl_in =\t\t %.2f \tmm"%(l_in))

print("\nl_n = \t\t%.2f \tmm"%(l_n))

print("\nl_s = \t\t%.2f \tmm"%(l_s))

print("\nR_s = \t\t%.2f \tmm"%(R_s))

print("\nRe_in = \t%d"%(Re_in))

m_ideal = np.pi*R_n**2*1e-6*Mi_def*np.sqrt(2*Rho_p*DeltaP*1e5)*1e3

print("\nm_ideal = \t%.2f \tg/s  (%.2f %%)"%(m_ideal,m_ideal/m_target*100))



print("\nInjetor Real:")

print("\nLambda = \t%.4f"%(Lambda))

print("\nA_eq = \t\t%.2f"%(A_eq))
        
print("\nMi_eq = \t%.4f"%(Mi_eq))

print("\nAlpha2_eq = \t%.1f \tdeg"%(2*Alpha2_eq*180/np.pi))

print("\nAlpha_in = \t%.1f \tdeg"%(Alpha_in))

print("\nKsi_in = \t%.4f"%(Ksi_in))

print("\nKsi_i = \t%.4f"%(Ksi_i))

print("\nMi_i = \t\t%.4f"%(Mi_i))

print("\nR_n1 = \t\t%.2f \tmm"%(R_n1))

m_real = np.pi*R_n**2*1e-6*Mi_i*np.sqrt(2*Rho_p*DeltaP*1e5)*1e3

print("\nm_real = \t%.2f \tg/s  (%.2f %%)"%(m_real,m_real/m_target*100))

print("\nNão Uniformidade Esperada = \t%.2f %%"%(I))




print("\nlen E: ",len(E))
print("len it: ",len(np.arange(1,Iteration)))

'''
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
'''

print("\n\n\n")





#####     #####     #####     #####     #####     #####     #####     #####
#####     #####     #####     Estágio de Combustível        #####     #####
#####     #####     #####     #####     #####     #####     #####     #####





#####  Fluid Properties: Start  #####

Propellant = "Ethanol 95%"



if(Propellant == "LOx"):

    Rho_p = 1141.0 #kg/m³
    
    Din_visc_p = 2.21e-3 #Pa.s
    
    m_target = 40 #g/s


elif(Propellant == "Ethanol 95%"):


    Rho_p = 800.0 #kg/m³
    
    Din_visc_p = 0.001095 #Pa.s

    m_target = 25 #g/s

elif(Propellant == "H2O 20ºC"):
    
    Rho_p = 965.31  #kg/m³

    Din_visc_p = 0.0008921 #Pa.s
    
    m_target = 50 #g/s



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

#fig, ax1 = plt.subplots(figsize=(14, 7))

    
for a in A:

    Mi = []
    for phi in Phi:
        
        
        Mi.append(1/(np.sqrt(a**2/(1-phi)+1/phi**2)))
    
    
    xy = ( Phi[Mi.index(max(Mi))]+0.03 , max(Mi)+0.004 )
    
    #plt.plot(Phi,Mi,color = 'black')
    #ax1.annotate("A = %.1f"%(a),xy)


Mi_mf = []


for phi in Phi:
    
    Mi_mf.append(phi*np.sqrt(phi/(2-phi)))

'''
plt.plot(Phi,Mi_mf, '-.', color = 'black')

plt.xticks(np.arange(0.0,1.0+0.1,0.05))
plt.xlim(0.0,1.0)

plt.yticks(np.arange(0.0,0.5+0.15,0.05))
plt.ylim(0.0,0.65)

plt.ylabel("μ")
plt.xlabel("φ")

plt.grid()

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


Alpha2_mf = []

for phi in Phi_mf:
    
    Alpha2_mf.append(2*180/np.pi*np.arctan(np.sqrt(2*(1-phi)/(phi))))
    


  
fig, ax2 = plt.subplots(figsize=(14, 7))

plt.plot(A,Alpha2_mf, color = 'black')
'''
plt.xticks(np.arange(0.0,25.0+2.0,1.0))
plt.xlim(0.0,25.00)

plt.yticks(np.arange(0.0,0.8+0.15,0.05))
plt.ylim(0.0,0.65)

plt.ylabel("μ")
plt.xlabel("A")
'''
plt.grid()

plt.show()



'''
fig, ax3 = plt.subplots(figsize=(14, 7))

plt.plot(A,Phi_mf, color = 'black')

plt.xticks(np.arange(0.0,A_max+0.5,1.0))
plt.xlim(0.0,A_max)

plt.yticks(np.arange(0.0,1.0+0.05,0.05))
plt.ylim(0.0,1.0)

plt.ylabel("φ")
plt.xlabel("A")

plt.grid()

plt.show()
'''


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


DeltaP = 5 #bar


Alpha2 = 110 #deg


for i in range(len(A)):
    
    if(abs(Alpha2_mf[i] - Alpha2) < 1.0):
        A_def = A[i]
        
print("A_def = ",A_def)


C = 2.00


N_in = 4


E = []



    
R_n = -2

R_n1 = -1

R_in = -1

r_in = -1

Iteration = 1

It_max = 1000

for i in range(len(A)):

    if(abs(A_def-A[i])<0.01):

        Mi_def = Mi_mf[i]

R_n =  0.475*np.sqrt(m_target*1e-3/(Mi_def*np.sqrt(Rho_p*DeltaP*1e5)))*1000 #mm


for i in range(len(A)):

    if(abs(A_def-A[i])<0.01):

        Phi_def = Phi_mf[i]





while(abs(R_n-R_n1)>1e-8 and Iteration <= It_max):

    
    
    if(R_n1 > 0):
        R_n = R_n1
        
        C = A_def*N_in*r_in**2/R_n
        
        R_in = R_n * C
        
        A_def = R_in*R_n/(N_in*r_in**2)
    
    
    Alpha2 = np.arctan(np.sqrt(2*(1-Phi_eq)/(Phi_eq)))
    
    
    

    
    r_in = np.sqrt(R_in*R_n/(N_in*A_def))
    
    
    
    l_in = 4*r_in #3~6
    
    
    
    l_n = 2*R_n #0.5~2
    
    
    
    l_s = 3*R_in #>2
    
    
    
    R_s = R_in+r_in
    
    
    
    Re_in = (2/np.pi)*m_target*1e-3/(N_in*r_in*1e-3*Din_visc_p)
    
    
    
    Lambda = 0.3164/(Re_in**0.25)
    
    
    
    ##### Recalc #####
    
    A_eq = R_in*R_n/(N_in*r_in**2+Lambda/2*R_in*(R_in-R_n))
    
    
    
    for i in range(len(A)):
    
        if(abs(A_eq-A[i])<0.01):
    
            Mi_eq = Mi_mf[i]
            
            
    
    
    for i in range(len(A)):
    
        if(abs(A_eq-A[i])<0.01):
    
            Phi_eq = Phi_mf[i]
    
    Alpha2_eq = np.arctan(np.sqrt(2*(1-Phi_eq)/(Phi_eq)))
    
    
    
    Alpha_in = 90 - 180/np.pi * np.arctan(R_s/l_in)
    
    
    
    Ksi_in = -1*Alpha_in/150 + 1.1
    
    
    
    Ksi_i = Ksi_in + Lambda*l_in/(2*r_in)
    
    
    
    Mi_i = Mi_eq/(np.sqrt(1+Ksi_i*Mi_eq**2*a**2/C**2))
    
    
    
    R_n1 = 0.475*np.sqrt(m_target*1e-3/(Mi_i*np.sqrt(Rho_p*DeltaP*1e5)))*1000 #mm
    
    
    I = 23.7/((R_in/R_n)**2.7*N_in**1.34*Phi_eq**1.1*(l_s/(2*R_s))**0.15)

    
    E.append(abs(r_in))#E.append(abs(R_n-R_n1))
    
    Iteration += 1
        
        
     
        
        



print('\033[4m' + '\t Estágio de Combustível \t\n\n' + '\033[0m')
        

print("Iteração: ",Iteration-1)

print("E_i = ",abs(R_n-R_n1))

print("\nInjetor ideal:")

print("\nA = \t\t%.2f"%(A_def))

print("\nMi = \t\t%.4f"%(Mi_def))

print("\nAlpha2 = \t%.1f \tdeg"%(2*Alpha2*180/np.pi))

print("\nR_n = \t\t%.2f \tmm"%(R_n))

print("\nD_n = \t\t%.2f \tmm"%(2*R_n))

print("\nC = \t\t%.2f"%(C))

print("\nR_in = \t\t%.2f \tmm"%(R_in))

print("\nN_in = \t\t%d"%(N_in))

print("\nr_in = \t\t%.2f \tmm"%(r_in))

print("\nd_in = \t\t%.2f \tmm"%(2*r_in))

print("\nl_in =\t\t %.2f \tmm"%(l_in))

print("\nl_n = \t\t%.2f \tmm"%(l_n))

print("\nl_s = \t\t%.2f \tmm"%(l_s))

print("\nR_s = \t\t%.2f \tmm"%(R_s))

print("\nRe_in = \t%d"%(Re_in))

m_ideal = np.pi*R_n**2*1e-6*Mi_def*np.sqrt(2*Rho_p*DeltaP*1e5)*1e3

print("\nm_ideal = \t%.2f \tg/s  (%.2f %%)"%(m_ideal,m_ideal/m_target*100))


print("\nInjetor Real:")

print("\nLambda = \t%.4f"%(Lambda))

print("\nA_eq = \t\t%.2f"%(A_eq))
        
print("\nMi_eq = \t%.4f"%(Mi_eq))

print("\nAlpha2_eq = \t%.1f \tdeg"%(2*Alpha2_eq*180/np.pi))

print("\nAlpha_in = \t%.1f \tdeg"%(Alpha_in))

print("\nKsi_in = \t%.4f"%(Ksi_in))

print("\nKsi_i = \t%.4f"%(Ksi_i))

print("\nMi_i = \t\t%.4f"%(Mi_i))

print("\nR_n1 = \t\t%.2f \tmm"%(R_n1))

m_real = np.pi*R_n**2*1e-6*Mi_i*np.sqrt(2*Rho_p*DeltaP*1e5)*1e3

print("\nm_real = \t%.2f \tg/s  (%.2f %%)"%(m_real,m_real/m_target*100))

print("\nNão Uniformidade Esperada = \t%.2f %%"%(I))






print("\nlen E: ",len(E))
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














