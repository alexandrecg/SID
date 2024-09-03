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


DeltaP_Ox = 5 #bar


Alpha2_Ox = 120 #deg

for i in range(len(A)):
    
    if(abs(Alpha2_mf[i] - Alpha2_Ox) < 0.5):
        A_Ox = A[i]
        
print("A_Ox = ",A_Ox)


##### ##### #####

C_Ox = 2.00


N_in_Ox = 4


E_Ox = []

    
R_n_Ox = -2

R_n1_Ox = -1

R_in_Ox = -1

r_in_Ox = -1

Iteration_Ox = 1

It_max_Ox = 10

for i in range(len(A)):

    if(abs(A_Ox-A[i])<0.01):

        Mi_Ox = Mi_mf[i]

R_n_Ox =  0.475*np.sqrt(m_target_Ox*1e-3/(Mi_Ox*np.sqrt(Rho_Ox*DeltaP_Ox*1e5)))*1000 #mm
            

##### ##### #####

deltaP_F_0 = 5 #bar


Alpha2_F_0 = 110 #deg


for i in range(len(A)):
    
    if(abs(Alpha2_mf[i] - Alpha2_F_0) < 1.0):
        A_F = A[i]
        
print("A_F = ",A_F)


C_F = 1.20


N_in_F = 6


E_F = []



Iteration_F = 1

It_max_F = 10

##### ##### #####



for i in range(len(A)):

    if(abs(A_F-A[i])<0.01):

        Mi_F = Mi_mf[i]


for i in range(len(A)):

    if(abs(A_F-A[i])<0.01):

        Phi_F = Phi_mf[i]


deltaP_F = -1
deltaP1_F = -2




delta_w = 1.00 #mm

delta_R = 0.5 #mm >= 0.3 mm

#Valores pre-definidos

R_n1_F = -1

R_in_F = -1

r_in_F = -1




# Main loop Start
while(abs(deltaP_F-deltaP_F_0)>1e-8 and Iteration_F <= It_max_F):
    
    print("It_Ox: ", Iteration_Ox)
    
    if(deltaP1_F < deltaP_F_0):
        C_Ox += 0.10
        
    elif(deltaP1_F < deltaP_F_0):
        C_Ox -= 0.10
    
    Iteration_Ox = 0
    
    #####     #####     #####     #####     #####     #####     #####     #####
    #####     #####     #####     Estágio de Oxidante           #####     #####
    #####     #####     #####     #####     #####     #####     #####     #####
    
    # Secondary loop Start - Oxidizer Stage
    while(abs(R_n_Ox-R_n1_Ox)>1e-8 and Iteration_Ox <= It_max_Ox):
        
        R_in_Ox = R_n_Ox * C_Ox
        
        if(R_n1_Ox > 0):
            R_n_Ox = R_n1_Ox
            A_Ox = R_in_Ox*R_n_Ox/(N_in_Ox*r_in_Ox**2)
        
        
        
        r_in_Ox = np.sqrt(R_in_Ox*R_n_Ox/(N_in_Ox*A_Ox))
        
        
        
        l_in_Ox = 4*r_in_Ox #3~6
        
        
        
        l_n_Ox = 2*R_n_Ox #0.5~2
        
        
        
        l_s_Ox = 3*R_in_Ox #>2
        
        
        
        R_s_Ox = R_in_Ox+r_in_Ox
        
        
        
        Re_in_Ox = (2/np.pi)*m_target_Ox*1e-3/(N_in_Ox*r_in_Ox*1e-3*Din_visc_Ox)
        
        
        
        Lambda_Ox = 0.3164/(Re_in_Ox**0.25)
        
        
        
        ##### Real #####
        
        A_eq_Ox = R_in_Ox*R_n_Ox/(N_in_Ox*r_in_Ox**2+Lambda_Ox/2*R_in_Ox*(R_in_Ox-R_n_Ox))
        
        
        
        for i in range(len(A)):
        
            if(abs(A_Ox-A[i])<0.01):
        
                Mi_eq_Ox = Mi_mf[i]
                
                
        
        
        for i in range(len(A)):
        
            if(abs(A_Ox-A[i])<0.01):
        
                Phi_eq_Ox = Phi_mf[i]
        
        Alpha2_eq_Ox = np.arctan(np.sqrt(2*(1-Phi_eq_Ox)/(Phi_eq_Ox)))
        
        
        
        Alpha_in_Ox = 90 - 180/np.pi * np.arctan(R_s_Ox/l_in_Ox)
        
        
        
        Ksi_in_Ox = -1*Alpha_in_Ox/150 + 1.1
        
        
        
        Ksi_i_Ox = Ksi_in_Ox + Lambda_Ox*l_in_Ox/(2*r_in_Ox)
        
        
        
        Mi_i_Ox = Mi_eq_Ox/(np.sqrt(1+Ksi_i_Ox*Mi_eq_Ox**2*A_Ox**2/C_Ox**2))
        
        
        
        R_n1_Ox = 0.475*np.sqrt(m_target_Ox*1e-3/(Mi_i_Ox*np.sqrt(Rho_Ox*DeltaP_Ox*1e5)))*1000 #mm
        
        I_Ox = 23.7/((R_in_Ox/R_n_Ox)**2.7*N_in_Ox**1.34*Phi_eq_Ox**1.1*(l_s_Ox/(2*R_s_Ox))**0.15)
        
        
        E_Ox.append(abs(A_Ox))
        
        Iteration_Ox += 1
    
    
 
    #####     #####     #####     #####     #####     #####     #####     #####
    #####     #####     #####     Interface dos Estágios        #####     #####
    #####     #####     #####     #####     #####     #####     #####     #####
    
    

    
    R_n_F = R_n_Ox + delta_w + delta_R
    

    #####     #####     #####     #####     #####     #####     #####     #####
    #####     #####     #####     Estágio de Combustível        #####     #####
    #####     #####     #####     #####     #####     #####     #####     #####
    

    

    if(Iteration_F > 1):
        deltaP_F = deltaP1_F
        
        C_F = A_F*N_in_F*r_in_F**2/R_n_F**2
        
    
    
    Alpha2_F = np.arctan(np.sqrt(2*(1-Phi_F)/(Phi_F)))
    
    
    R_in_F = R_n_F * C_F

    
    r_in_F = np.sqrt(R_in_F*R_n_F/(N_in_F*A_F))
    
    
    
    l_in_F = 4*r_in_F #3~6
    
    
    
    l_n_F = 2*R_n_F #0.5~2
    
    
    
    l_s_F = 3*R_in_F #>2
    
    
    
    R_s_F = R_in_F+r_in_F
    
    
    
    Re_in_F = (2/np.pi)*m_target_F*1e-3/(N_in_F*r_in_F*1e-3*Din_visc_F)
    
    deltaP_F = 0.05*m_target_F**2/(Mi_F**2*Rho_F*R_n_F**4)
    
    Lambda_F = 0.3164/(Re_in_F**0.25)
    
    
    
    ##### Real #####
    
    A_eq_F = R_in_F*R_n_F/(N_in_F*r_in_F**2+Lambda_F/2*R_in_F*(R_in_F-R_n_F))
    
    
    
    for i in range(len(A)):
    
        if(abs(A_eq_F-A[i])<0.01):
    
            Mi_eq_F = Mi_mf[i]
            
            
    
    
    for i in range(len(A)):
    
        if(abs(A_eq_F-A[i])<0.01):
    
            Phi_eq_F = Phi_mf[i]
    
    Alpha2_eq_F = np.arctan(np.sqrt(2*(1-Phi_eq_F)/(Phi_eq_F)))
    
    
    
    Alpha_in_F = 90 - 180/np.pi * np.arctan(R_s_F/l_in_F)
    
    
    
    Ksi_in_F = -1*Alpha_in_F/150 + 1.1
    
    
    
    Ksi_i_F = Ksi_in_F + Lambda_F*l_in_F/(2*r_in_F)
    
    
    
    Mi_i_F = Mi_eq_F/(np.sqrt(1+Ksi_i_F*Mi_eq_F**2*A_F**2/C_F**2))
    
    
    
    deltaP1_F = 1/(2*np.pi**2)*(m_target_F*1e-3)**2/(Mi_i_F**2*Rho_F*(R_n_F*1e-3)**4)*1e-5 #bar

    print("\ndeltaP1_F = ", deltaP1_F, " bar")
    
    I_F = 23.7/((R_in_F/R_n_F)**2.7*N_in_F**1.34*Phi_eq_F**1.1*(l_s_F/(2*R_s_F))**0.15)

    
    E_F.append(abs(C_F))#E.append(abs(R_n-R_n1))
    
    
    
    
    
    Iteration_F += 1



    #####     #####     #####     #####     #####     #####     #####     #####
    #####     #####     #####     Relatório           #####     #####     #####  
    #####     #####     #####     #####     #####     #####     #####     #####     



print('\033[4m' + '\t Estágio de Oxidante \t\n\n' + '\033[0m')
        

#print("Iteração: ",Iteration_Ox-1)

print("E_i = ",abs(R_n_Ox-R_n1_Ox))



print("\nInjetor Ideal:")

print("\nA = \t\t%.2f"%(A_Ox))

print("\nMi = \t\t%.4f"%(Mi_Ox))

print("\nR_n = \t\t%.2f \tmm"%(R_n_Ox))

print("\nD_n = \t\t%.2f \tmm"%(2*R_n_Ox))

print("\nC = \t\t%.2f"%(C_Ox))

print("\nR_in = \t\t%.2f \tmm"%(R_in_Ox))

print("\nN_in = \t\t%d"%(N_in_Ox))

print("\nr_in = \t\t%.2f \tmm"%(r_in_Ox))

print("\nd_in = \t\t%.2f \tmm"%(2*r_in_Ox))

print("\nl_in =\t\t %.2f \tmm"%(l_in_Ox))

print("\nl_n = \t\t%.2f \tmm"%(l_n_Ox))

print("\nl_s = \t\t%.2f \tmm"%(l_s_Ox))

print("\nR_s = \t\t%.2f \tmm"%(R_s_Ox))

print("\nRe_in = \t%d"%(Re_in_Ox))

m_ideal_Ox = np.pi*R_n_Ox**2*1e-6*Mi_Ox*np.sqrt(2*Rho_Ox*DeltaP_Ox*1e5)*1e3

print("\nm_ideal = \t%.2f \tg/s  (%.2f %%)"%(m_ideal_Ox,m_ideal_Ox/m_target_Ox*100))



print("\nInjetor Real:")

print("\nLambda = \t%.4f"%(Lambda_Ox))

print("\nA_eq = \t\t%.2f"%(A_eq_Ox))
        
print("\nMi_eq = \t%.4f"%(Mi_eq_Ox))

print("\nAlpha2_eq = \t%.1f \tdeg"%(2*Alpha2_eq_Ox*180/np.pi))

print("\nAlpha_in = \t%.1f \tdeg"%(Alpha_in_Ox))

print("\nKsi_in = \t%.4f"%(Ksi_in_Ox))

print("\nKsi_i = \t%.4f"%(Ksi_i_Ox))

print("\nMi_i = \t\t%.4f"%(Mi_i_Ox))

print("\nR_n1 = \t\t%.2f \tmm"%(R_n1_Ox))

m_real_Ox = np.pi*R_n_Ox**2*1e-6*Mi_i_Ox*np.sqrt(2*Rho_Ox*DeltaP_Ox*1e5)*1e3

print("\nm_real = \t%.2f \tg/s  (%.2f %%)"%(m_real_Ox,m_real_Ox/m_target_Ox*100))

print("\nNão Uniformidade Esperada = \t%.2f %%"%(I_Ox))




print("\nlen E: ",len(E_Ox))
print("len it: ",len(np.arange(1,Iteration_Ox)))


print("\n\n\n")


print('\033[4m' + '\t Estágio de Combustível \t\n\n' + '\033[0m')
        

print("Iteração: ",Iteration_F-1)

print("E_i = ",abs(R_n_F-R_n1_F))

print("\nInjetor ideal:")

print("\nA = \t\t%.2f"%(A_F))

print("\nMi = \t\t%.4f"%(Mi_F))

print("\nAlpha2 = \t%.1f \tdeg"%(2*Alpha2_F*180/np.pi))

print("\nR_n = \t\t%.2f \tmm"%(R_n_F))

print("\nD_n = \t\t%.2f \tmm"%(2*R_n_F))

print("\nC = \t\t%.2f"%(C_F))

print("\nR_in = \t\t%.2f \tmm"%(R_in_F))

print("\nN_in = \t\t%d"%(N_in_F))

print("\nr_in = \t\t%.2f \tmm"%(r_in_F))

print("\nd_in = \t\t%.2f \tmm"%(2*r_in_F))

print("\nl_in =\t\t %.2f \tmm"%(l_in_F))

print("\nl_n = \t\t%.2f \tmm"%(l_n_F))

print("\nl_s = \t\t%.2f \tmm"%(l_s_F))

print("\nR_s = \t\t%.2f \tmm"%(R_s_F))

print("\nRe_in = \t%d"%(Re_in_F))

m_ideal_F = np.pi*R_n_F**2*1e-6*Mi_F*np.sqrt(2*Rho_F*deltaP_F*1e5)*1e3

print("\nm_ideal = \t%.2f \tg/s  (%.2f %%)"%(m_ideal_F,m_ideal_F/m_target_F*100))


print("\nInjetor Real:")

print("\nLambda = \t%.4f"%(Lambda_F))

print("\nA_eq = \t\t%.2f"%(A_eq_F))
        
print("\nMi_eq = \t%.4f"%(Mi_eq_F))

print("\nAlpha2_eq = \t%.1f \tdeg"%(2*Alpha2_eq_F*180/np.pi))

print("\nAlpha_in = \t%.1f \tdeg"%(Alpha_in_F))

print("\nKsi_in = \t%.4f"%(Ksi_in_F))

print("\nKsi_i = \t%.4f"%(Ksi_i_F))

print("\nMi_i = \t\t%.4f"%(Mi_i_F))

print("\nR_n1 = \t\t%.2f \tmm"%(R_n1_F))

m_real_F = np.pi*R_n_F**2*1e-6*Mi_i_F*np.sqrt(2*Rho_F*deltaP_F*1e5)*1e3

print("\nm_real = \t%.2f \tg/s  (%.2f %%)"%(m_real_F,m_real_F/m_target_F*100))

print("\nNão Uniformidade Esperada = \t%.2f %%"%(I_F))





print("\nlen E_F: ",len(E_F))
print("len it_F: ",len(np.arange(1,Iteration_F)))


print("\nlen E_Ox: ",len(E_Ox))
print("len it_Ox: ",len(np.arange(1,len(E_Ox)+1)))

fig, ax3 = plt.subplots(figsize=(14, 7))

plt.plot(np.arange(1,len(E_Ox)+1),E_Ox, color = 'black')

#plt.xticks(np.arange(0.0,A_max+0.5,1.0))
#plt.xlim(0.0,A_max)

#plt.yticks(np.arange(0.0,1.0+0.05,0.05))
#plt.ylim(0.0,1.0)

plt.ylabel("E_Ox")
plt.xlabel("it.")

plt.grid()



fig, ax4 = plt.subplots(figsize=(14, 7))

plt.plot(np.arange(1,Iteration_F),E_F, color = 'black')

#plt.xticks(np.arange(0.0,A_max+0.5,1.0))
#plt.xlim(0.0,A_max)

#plt.yticks(np.arange(0.0,1.0+0.05,0.05))
#plt.ylim(0.0,1.0)

plt.ylabel("E_F")
plt.xlabel("it.")

plt.grid()

plt.show()



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####




    

        
        
     
        
        


















