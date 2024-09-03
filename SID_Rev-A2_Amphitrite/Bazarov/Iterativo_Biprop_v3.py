# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 22:11:09 2020

@author: A. Goulart
"""

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt





#####  Fluid Properties: Start  #####

Oxidizer = "LOx"

Rho_Ox = 1141.0 #kg/m³

Din_visc_Ox = 2.21e-3 #Pa.s

m_target_Ox = 40 #g/s



Fuel = "Ethanol 95%"

Rho_F = 800.0 #kg/m³

Din_visc_F = 0.001095 #Pa.s

m_target_F = 25 #g/s


#Estágio Interno

Rho_Int = Rho_F

Din_visc_Int = Din_visc_F

m_target_Int = m_target_F

#Estágio Externo

Rho_Ext = Rho_Ox

Din_visc_Ext = Din_visc_Ox

m_target_Ext = m_target_Ox


#####  Fluid Properties: End  #####

Phi_min = 0.0
    
Phi_max = 1.0

Phi_step = (Phi_max-Phi_min)/1000

Phi = np.arange(Phi_min,Phi_max,Phi_step)

##### ##### #####

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
    


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


DeltaP_Int = 5 #bar


Alpha2_Int = 120 #deg

for i in range(len(A)):
    
    if(abs(Alpha2_mf[i] - Alpha2_Int) < 0.5):
        A_Int = A[i]
        
print("A_Int = ",A_Int)


##### ##### #####

C_Int = 2.00


N_in_Int = 4


E_Int = []

    
R_n_Int = -2

R_n1_Int = -1

R_in_Int = -1

r_in_Int = -1

Iteration_Int = 1

It_max_Int = 50

for i in range(len(A)):

    if(abs(A_Int-A[i])<0.01):

        Mi_Int = Mi_mf[i]


            

##### ##### #####

deltaP_Ext_0 = 5 #bar


Alpha2_Ext_0 = 110 #deg


for i in range(len(A)):
    
    if(abs(Alpha2_mf[i] - Alpha2_Ext_0) < 1.0):
        A_Ext = A[i]
        
print("A_Ext = ",A_Ext)


C_Ext = 1.20


N_in_Ext = 6


E_Ext = []



Iteration_Ext = 1

It_max_Ext = 50

##### ##### #####



for i in range(len(A)):

    if(abs(A_Ext-A[i])<0.01):

        Mi_Ext = Mi_mf[i]


for i in range(len(A)):

    if(abs(A_Ext-A[i])<0.01):

        Phi_Ext = Phi_mf[i]


deltaP_Ext = -1
deltaP1_Ext = -2




delta_w = 1.00 #mm

delta_R = 0.5 #mm >= 0.3 mm

#Valores pre-definidos

R_n1_Ext = -1

R_in_Ext = -1

r_in_Ext = -1


It_Int = []
It_Ext = []

Cint = []


# Main loop Start
while(abs(deltaP_Ext-deltaP_Ext_0)>1e-8 and Iteration_Ext <= It_max_Ext):
    
    print("It_Int: ", Iteration_Int)
    
    if(deltaP1_Ext < deltaP_Ext_0):
        C_Int -= 0.10
        
    elif(deltaP1_Ext < deltaP_Ext_0):
        C_Int += 0.10
    
    if(C_Int < 1.0): C_Int = 1.0
    
    Iteration_Int = 0
    
    #####     #####     #####     #####     #####     #####     #####     #####
    #####     #####     #####     Estágio de Oxidante           #####     #####
    #####     #####     #####     #####     #####     #####     #####     #####
    
    # Secondary loop Start - Internal Stage
    while(abs(R_n_Int-R_n1_Int)>1e-8 and Iteration_Int <= It_max_Int):
        
        
        
        if(R_n1_Int > 0):
            R_n_Int = R_n1_Int
            
            r_in_Int = np.sqrt(R_in_Int*R_n_Int/(N_in_Int*A_Int))
            
            A_Int = R_in_Int*R_n_Int/(N_in_Int*r_in_Int**2)
            
        else:
            R_n_Int =  0.475*np.sqrt(m_target_Int*1e-3/(Mi_Int*np.sqrt(Rho_Int*DeltaP_Int*1e5)))*1000 #mm
            
            r_in_Int = np.sqrt(R_in_Int*R_n_Int/(N_in_Int*A_Int))
        
        
        
        R_in_Int = R_n_Int * C_Int
        
        
        l_in_Int = 4*r_in_Int #3~6
        
        
        
        l_n_Int = 2*R_n_Int #0.5~2
        
        
        
        l_s_Int = 3*R_in_Int #>2
        
        
        
        R_s_Int = R_in_Int+r_in_Int
        
        
        
        Re_in_Int = (2/np.pi)*m_target_Int*1e-3/(N_in_Int*r_in_Int*1e-3*Din_visc_Int)
        
        
        
        Lambda_Int = 0.3164/(Re_in_Int**0.25)
        
        
        
        ##### Real #####
        
        A_eq_Int = R_in_Int*R_n_Int/(N_in_Int*r_in_Int**2+Lambda_Int/2*R_in_Int*(R_in_Int-R_n_Int))
        
        
        
        for i in range(len(A)):
        
            if(abs(A_Int-A[i])<0.01):
        
                Mi_eq_Int = Mi_mf[i]
                
                
        
        
        for i in range(len(A)):
        
            if(abs(A_Int-A[i])<0.01):
        
                Phi_eq_Int = Phi_mf[i]
        
        Alpha2_eq_Int = np.arctan(np.sqrt(2*(1-Phi_eq_Int)/(Phi_eq_Int)))
        
        
        
        Alpha_in_Int = 90 - 180/np.pi * np.arctan(R_s_Int/l_in_Int)
        
        
        
        Ksi_in_Int = -1*Alpha_in_Int/150 + 1.1
        
        
        
        Ksi_i_Int = Ksi_in_Int + Lambda_Int*l_in_Int/(2*r_in_Int)
        
        
        
        Mi_i_Int = Mi_eq_Int/(np.sqrt(1+Ksi_i_Int*Mi_eq_Int**2*A_Int**2/C_Int**2))
        
        
        
        R_n1_Int = 0.475*np.sqrt(m_target_Int*1e-3/(Mi_i_Int*np.sqrt(Rho_Int*DeltaP_Int*1e5)))*1000 #mm
        
        I_Int = 23.7/((R_in_Int/R_n_Int)**2.7*N_in_Int**1.34*Phi_eq_Int**1.1*(l_s_Int/(2*R_s_Int))**0.15)
        
        
        E_Int.append(abs(R_n_Int-R_n1_Int))
        
        Cint.append(C_Int)
        
        It_Int.append(Iteration_Int)
        It_Ext.append(Iteration_Ext)
        
        Iteration_Int += 1
    
    
 
    #####     #####     #####     #####     #####     #####     #####     #####
    #####     #####     #####     Interface dos Estágios        #####     #####
    #####     #####     #####     #####     #####     #####     #####     #####
    
    

    
    R_n_Ext = R_n_Int + delta_w + delta_R
    

    #####     #####     #####     #####     #####     #####     #####     #####
    #####     #####     #####     Estágio de Combustível        #####     #####
    #####     #####     #####     #####     #####     #####     #####     #####
    

    

    if(Iteration_Ext > 1):
        deltaP_Ext = deltaP1_Ext
        
        C_Ext = A_Ext*N_in_Ext*r_in_Ext**2/R_n_Ext**2
        
        if(C_Ext < 1.0): C_Ext = 1.0
        
    
    
    Alpha2_Ext = np.arctan(np.sqrt(2*(1-Phi_Ext)/(Phi_Ext)))
    
    
    R_in_Ext = R_n_Ext * C_Ext

    
    r_in_Ext = np.sqrt(R_in_Ext*R_n_Ext/(N_in_Ext*A_Ext))
    
    
    
    l_in_Ext = 4*r_in_Ext #3~6
    
    
    
    l_n_Ext = 2*R_n_Ext #0.5~2
    
    
    
    l_s_Ext = 3*R_in_Ext #>2
    
    
    
    R_s_Ext = R_in_Ext+r_in_Ext
    
    
    
    Re_in_Ext = (2/np.pi)*m_target_Ext*1e-3/(N_in_Ext*r_in_Ext*1e-3*Din_visc_Ext)
    
    deltaP_Ext = 0.05*m_target_Ext**2/(Mi_Ext**2*Rho_Ext*R_n_Ext**4)
    
    Lambda_Ext = 0.3164/(Re_in_Ext**0.25)
    
    
    
    ##### Real #####
    
    A_eq_Ext = R_in_Ext*R_n_Ext/(N_in_Ext*r_in_Ext**2+Lambda_Ext/2*R_in_Ext*(R_in_Ext-R_n_Ext))
    
    
    
    for i in range(len(A)):
    
        if(abs(A_eq_Ext-A[i])<0.01):
    
            Mi_eq_Ext = Mi_mf[i]
            
            
    
    
    for i in range(len(A)):
    
        if(abs(A_eq_Ext-A[i])<0.01):
    
            Phi_eq_Ext = Phi_mf[i]
    
    Alpha2_eq_Ext = np.arctan(np.sqrt(2*(1-Phi_eq_Ext)/(Phi_eq_Ext)))
    
    
    
    Alpha_in_Ext = 90 - 180/np.pi * np.arctan(R_s_Ext/l_in_Ext)
    
    
    
    Ksi_in_Ext = -1*Alpha_in_Ext/150 + 1.1
    
    
    
    Ksi_i_Ext = Ksi_in_Ext + Lambda_Ext*l_in_Ext/(2*r_in_Ext)
    
    
    
    Mi_i_Ext = Mi_eq_Ext/(np.sqrt(1+Ksi_i_Ext*Mi_eq_Ext**2*A_Ext**2/C_Ext**2))
    
    
    
    deltaP1_Ext = 1/(2*np.pi**2)*(m_target_Ext*1e-3)**2/(Mi_i_Ext**2*Rho_Ext*(R_n_Ext*1e-3)**4)*1e-5 #bar

    print("\ndeltaP1_Ext = ", deltaP1_Ext, " bar")
    
    I_Ext = 23.7/((R_in_Ext/R_n_Ext)**2.7*N_in_Ext**1.34*Phi_eq_Ext**1.1*(l_s_Ext/(2*R_s_Ext))**0.15)

    
    E_Ext.append(abs(C_Ext))#E.append(abs(R_n-R_n1))
    
    
    
    
    
    Iteration_Ext += 1



    #####     #####     #####     #####     #####     #####     #####     #####
    #####     #####     #####     Relatório           #####     #####     #####  
    #####     #####     #####     #####     #####     #####     #####     #####     



print('033[4m' + '\t Estágio Interno \t\n\n' + '033[0m')
        

#print("Iteração: ",Iteration_Int-1)

print("E_i = ",abs(R_n_Int-R_n1_Int))



print("\nInjetor Ideal:")

print("\nA = \t\t%.2f"%(A_Int))

print("\nMi = \t\t%.4f"%(Mi_Int))

print("\nR_n = \t\t%.2f \tmm"%(R_n_Int))

print("\nD_n = \t\t%.2f \tmm"%(2*R_n_Int))

print("\nC = \t\t%.2f"%(C_Int))

print("\nR_in = \t\t%.2f \tmm"%(R_in_Int))

print("\nN_in = \t\t%d"%(N_in_Int))

print("\nr_in = \t\t%.2f \tmm"%(r_in_Int))

print("\nd_in = \t\t%.2f \tmm"%(2*r_in_Int))

print("\nl_in =\t\t %.2f \tmm"%(l_in_Int))

print("\nl_n = \t\t%.2f \tmm"%(l_n_Int))

print("\nl_s = \t\t%.2f \tmm"%(l_s_Int))

print("\nR_s = \t\t%.2f \tmm"%(R_s_Int))

print("\nRe_in = \t%f"%(Re_in_Int))

m_ideal_Int = np.pi*R_n_Int**2*1e-6*Mi_Int*np.sqrt(2*Rho_Int*DeltaP_Int*1e5)*1e3

print("\nm_ideal = \t%.2f \tg/s  (%.2f %%)"%(m_ideal_Int,m_ideal_Int/m_target_Int*100))



print("\nInjetor Real:")

print("\nLambda = \t%.4f"%(Lambda_Int))

print("\nA_eq = \t\t%.2f"%(A_eq_Int))
        
print("\nMi_eq = \t%.4f"%(Mi_eq_Int))

print("\nAlpha2_eq = \t%.1f \tdeg"%(2*Alpha2_eq_Int*180/np.pi))

print("\nAlpha_in = \t%.1f \tdeg"%(Alpha_in_Int))

print("\nKsi_in = \t%.4f"%(Ksi_in_Int))

print("\nKsi_i = \t%.4f"%(Ksi_i_Int))

print("\nMi_i = \t\t%.4f"%(Mi_i_Int))

print("\nR_n1 = \t\t%.2f \tmm"%(R_n1_Int))

m_real_Int = np.pi*R_n_Int**2*1e-6*Mi_i_Int*np.sqrt(2*Rho_Int*DeltaP_Int*1e5)*1e3

print("\nm_real = \t%.2f \tg/s  (%.2f %%)"%(m_real_Int,m_real_Int/m_target_Int*100))

print("\nNão Uniformidade Esperada = \t%.2f %%"%(I_Int))




print("\nlen E: ",len(E_Int))
print("len it: ",len(np.arange(1,Iteration_Int)))


print("\n\n\n")


print('033[4m' + '\t Estágio Externo \t\n\n' + '033[0m')
        

print("Iteração: ",Iteration_Ext-1)

print("E_i = ",abs(R_n_Ext-R_n1_Ext))

print("\nInjetor ideal:")

print("\nA = \t\t%.2f"%(A_Ext))

print("\nMi = \t\t%.4f"%(Mi_Ext))

print("\nAlpha2 = \t%.1f \tdeg"%(2*Alpha2_Ext*180/np.pi))

print("\nR_n = \t\t%.2f \tmm"%(R_n_Ext))

print("\nD_n = \t\t%.2f \tmm"%(2*R_n_Ext))

print("\nC = \t\t%.2f"%(C_Ext))

print("\nR_in = \t\t%.2f \tmm"%(R_in_Ext))

print("\nN_in = \t\t%d"%(N_in_Ext))

print("\nr_in = \t\t%.2f \tmm"%(r_in_Ext))

print("\nd_in = \t\t%.2f \tmm"%(2*r_in_Ext))

print("\nl_in =\t\t %.2f \tmm"%(l_in_Ext))

print("\nl_n = \t\t%.2f \tmm"%(l_n_Ext))

print("\nl_s = \t\t%.2f \tmm"%(l_s_Ext))

print("\nR_s = \t\t%.2f \tmm"%(R_s_Ext))

print("\nRe_in = \t%d"%(Re_in_Ext))

m_ideal_Ext = np.pi*R_n_Ext**2*1e-6*Mi_Ext*np.sqrt(2*Rho_Ext*deltaP_Ext*1e5)*1e3

print("\nm_ideal = \t%.2f \tg/s  (%.2f %%)"%(m_ideal_Ext,m_ideal_Ext/m_target_Ext*100))


print("\nInjetor Real:")

print("\nLambda = \t%.4f"%(Lambda_Ext))

print("\nA_eq = \t\t%.2f"%(A_eq_Ext))
        
print("\nMi_eq = \t%.4f"%(Mi_eq_Ext))

print("\nAlpha2_eq = \t%.1f \tdeg"%(2*Alpha2_eq_Ext*180/np.pi))

print("\nAlpha_in = \t%.1f \tdeg"%(Alpha_in_Ext))

print("\nKsi_in = \t%.4f"%(Ksi_in_Ext))

print("\nKsi_i = \t%.4f"%(Ksi_i_Ext))

print("\nMi_i = \t\t%.4f"%(Mi_i_Ext))

print("\nR_n1 = \t\t%.2f \tmm"%(R_n1_Ext))

m_real_Ext = np.pi*R_n_Ext**2*1e-6*Mi_i_Ext*np.sqrt(2*Rho_Ext*deltaP_Ext*1e5)*1e3

print("\nm_real = \t%.2f \tg/s  (%.2f %%)"%(m_real_Ext,m_real_Ext/m_target_Ext*100))

print("\nNão Uniformidade Esperada = \t%.2f %%"%(I_Ext))





print("\nlen E_Ext: ",len(E_Ext))
print("len it_Ext: ",len(np.arange(1,Iteration_Ext)))


print("\nlen E_Int: ",len(E_Int))
print("len it_Int: ",len(np.arange(1,len(E_Int)+1)))


print("\nClen Cint: ", len(Cint))
print("len It_Int: ", len(It_Int))
print("len It_Ext: ", len(It_Ext))

fig, ax1 = plt.subplots(figsize=(14, 7))

plt.plot(np.arange(1,len(E_Int)+1),E_Int, color = 'black')

#plt.xticks(np.arange(0.0,A_max+0.5,1.0))
#plt.xlim(0.0,A_max)

#plt.yticks(np.arange(0.0,1.0+0.05,0.05))
#plt.ylim(0.0,1.0)

plt.ylabel("E_Int")
plt.xlabel("it.")

plt.grid()

########## ##########

fig, ax2 = plt.subplots(figsize=(14, 7))

plt.plot(np.arange(1,len(E_Int)+1),Cint, color = 'b')

#plt.xticks(np.arange(0.0,A_max+0.5,1.0))
#plt.xlim(0.0,A_max)

#plt.yticks(np.arange(0.0,1.0+0.05,0.05))
#plt.ylim(0.0,1.0)

plt.ylabel("C_Int")
plt.xlabel("it.")

plt.grid()

########## ##########

fig, ax3 = plt.subplots(figsize=(14, 7))

plt.plot(np.arange(1,Iteration_Ext),E_Ext, color = 'black')

#plt.xticks(np.arange(0.0,A_max+0.5,1.0))
#plt.xlim(0.0,A_max)

#plt.yticks(np.arange(0.0,1.0+0.05,0.05))
#plt.ylim(0.0,1.0)

plt.ylabel("E_Ext")
plt.xlabel("it.")

plt.grid()

########## ##########

fig, ax4 = plt.subplots(figsize=(14, 7))

plt.plot(It_Int,It_Ext, color = 'black')

#plt.xticks(np.arange(0.0,A_max+0.5,1.0))
#plt.xlim(0.0,A_max)

#plt.yticks(np.arange(0.0,1.0+0.05,0.05))
#plt.ylim(0.0,1.0)

plt.ylabel("It_Ext")
plt.xlabel("It_Int")

plt.grid()

plt.show()


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####




    

        
        
     
        
        


















