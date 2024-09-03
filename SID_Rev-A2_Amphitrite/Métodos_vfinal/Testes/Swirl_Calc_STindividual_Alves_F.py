# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 18:12:02 2019

@author: Alexandre Goulart
"""

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt


#SIFG - Swirl Injector from Geometry

#####  Fluid Properties  #####

Fluid = "F" # "F" ou "Ox"



if(Fluid == "F"):
    Propellant = "RP-1"
    
    Rho_p = 800 #800.0 #kg/m³
    
    Din_visc_p = 1.92e-3#1.447*1e-3 #0.001095 #Pa.s
    
    m_target = 72.2#64.8 #g/s


elif(Fluid == "Ox"):
    Propellant = "LOx"
    
    Rho_p = 1141.0 #kg/m³
    
    Din_visc_p = 2.21e-3 #Pa.s
    
    m_target = 172.9 #g/s






print("m_target: \t%.2f \tg/s\n"%(np.round(m_target,2)))

#####  Operation Parameters  #####

Delta_P0 = 0.710 #MPa


#####  Injector Geometry  #####

print("\n\nGeometria do injetor")

N_in = 6

print("Número de canais de entrada (N_in): \t\t%.1f"%(N_in))

d_in = 2*0.40 #mm

print("Diâmetro dos canais de entrada (d_in): \t\t%.2f \tmm"%(d_in))

r_in = d_in/2 #mm

#k_L = 2.0 #2.0 a 3.0

#Fillet = 1.1* d_in#0.20 #mm (Chamfro de entrada)

L_in = 2.4 #mm #k_L*d_in

print("Comprimento dos canais de entrada (L_in): \t%.2f \tmm"%(L_in))

D_cv = 2*5.10 #mm

print("Diâmetro da Câmara de Vórtice (D_cv): \t\t%.2f \tmm"%(D_cv))

R_cv = D_cv/2 #mm

L_cv = 10.5 #mm

print("Comprimento da Câmara de Vórtice (L_cv): \t%.2f \tmm"%(L_cv))

D_n = 2*5.10 #mm ###################################

print("Diâmetro do Orifício de saída (D_n): \t\t%.2f \tmm"%(D_n))

R_n = D_n/2 #mm

R_in = 4.7 #R_cv-r_in #mm ##################################

print("Raio de Entrada dos canais (R_in): \t\t%.2f \tmm"%(R_in))

A_n = np.pi*R_n**2 #mm²

Fillet = 0.91 #mm

##### ##### ##### funções estágio ideal ##### ##### #####

"""
def calc_A_i(R_in)
    A_i = R_in*R_n/(N_in*(r_in)**2)
 
"""

##### ##### ##### funções estágio real ##### ##### #####


print("\nDelta_P0: \t\t%.2f bar"%(Delta_P0*10))



print("\n\nDados do Estágio ideal")

### Ideal Injector Parameters  ###


#parâmetro geometrico

def calc_A (R_in, R_n, N_in, r_in, Lambda):
    return R_in*1e-3*R_n*1e-3/(N_in*(r_in*1e-3)**2+Lambda/2*R_in*1e-3*(R_in*1e-3-R_n*1e-3))

A_i = calc_A (R_in, R_n, N_in, r_in, 0) 

## Phi_i calculation ##
 # Start #


def calc_Phi(A):
        
    Phi_0 = 0.1 #initial value
    
    def f_Phi(x):
        return (1-x)*np.sqrt(2)/x**1.5-A
    
    Phi = opt.fsolve(f_Phi, Phi_0)
    
    if(len(Phi) == 1):
        Phi = Phi[0]
        return Phi
        
    elif(len(Phi) == 0):
        print("ERROR: no Phi value found")
        
    else:
        print("ERROR: Phi multiple values")

Phi_i = calc_Phi(A_i)

 # End #
 
#Coeficiente de Descarga ideal
Mu_i = Phi_i*np.sqrt(Phi_i/(2-Phi_i))

#Espessura de camada líquida ideal
t_as_i = R_n*(1-np.sqrt(1-Phi_i)) #mm

if(t_as_i <= 0):
    print("ERROR: invalid t_as value")
    
elif(t_as_i <= 0.1):
    print("Caution! Low value of t_as")
    
#Vazão mássica ideal
m_i = Mu_i*A_n*1e-6*np.sqrt(2*Delta_P0*1e6*Rho_p)*1e3 #vazão mássica ideal [g/s]


#Coeficiente de redução de área do estágio
C = R_in/R_n


print("Parâmetro Geométrico ideal (A_i): \t%.2f"%(A_i))
print("Coef. de área livre ideal (Phi_i): \t%.4f"%(Phi_i))
print("Coeficiente de descarga ideal(Mu_i): \t%.4f"%(Mu_i))
print("Espessura da camada líquida (t_as_i): \t%.4f \tmm"%(t_as_i))
print("Vazão mássica ideal (m_i): \t\t%.1f \tg/s"%(m_i))

##### Cálulo do ângulo de spray ideal #####

Alpha2 = 2*np.arctan(np.sqrt(2*(1-Phi_i)/Phi_i))*180/np.pi  #2*180/np.pi*np.arcsin(2*Mu_i*A_i/((1+np.sqrt(1-Phi_i))*np.sqrt(1-0*Mu_i**2*A_i**2/C**2)))

print("Ângulo de Abertura do Spray (Alpha2): \t%.1f \tdeg"%(Alpha2))

#####   Cálculo do Estágio Real   #####

print("\nDados do Estágio real")


#Velocidade Tangencial
W = m_i*1e-3/(N_in*np.pi*(r_in*1e-3)**2*Rho_p) #Tangential Velocity



#Número de Reynolds nos canais tangenciais
Re_in = W*d_in*1e-3/(Din_visc_p/Rho_p)



if(Re_in < 2600): print("ERROR: Low Reynolds number (%.2E < 2600)"%(Re_in))
elif(Re_in > 5e4): print("Caution: Reynolds number over 5e4")

#Coeficiente de resistência de superfície
Lambda = 0.3164*Re_in**(-0.25)



#Parâmetro Geométrico equivalente
A_eq = calc_A (R_in, R_n, N_in, r_in, Lambda) 



##### Cálculo do Coeficiente de área livre real

Phi_eq = calc_Phi(A_eq)


#####

#Coeficiente de Descarga
    
Mu_w = Phi_eq**1.5/(np.sqrt(2-Phi_eq))


#Espessura de camada líquida real
t_as_w = R_n*(1-np.sqrt(1-Phi_eq)) #mm



if(t_as_w <= 0):
    print("ERROR: invalid fluid layer thickness value")
    
elif(t_as_w <= 0.1):
    print("Caution! Low fluid layer thickness")


#Coeficiente de perda de momento angular
K = A_eq/A_i



##### Cálculo da perda nos canais tangenciais #####

Alpha_in = 90 - 180/np.pi*np.arctan(R_cv/L_in)

#Coeficiente de perda
Ksi_in = -0.0204*(Fillet/d_in)**4-0.026*(Fillet/d_in)**3+0.3587*(Fillet/d_in)**2-0.6995*(Fillet/d_in)+0.5002
#Ksi_in = -1*Alpha_in/150 + 1.1


#Fator de perda nos canais tangenciais
Ksi = Ksi_in + Lambda*L_in/(2*r_in)

 #coeficiente de descarga total
Mu_i = 1/np.sqrt(1/Phi_eq**2+A_i**2/(1-Phi_eq)+(Ksi*N_in)*A_i**2/C**2)
#Mu_i = Mu_w/(np.sqrt(1+Mu_w**2*Ksi*A_i**2/C**2))

#Vazão mássica real
m_r = Mu_i*np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P0*1e6)*1000 #np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P*1e6)/np.sqrt((2-Phi_eq)/Phi_eq**3+Ksi*A_i**2/R_in**2)*1e3
#m_r = Mu_i*np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P*1e6)*1e3

print("Velocidade tangencial (W): \t\t%.2f \tm/s"%(W))
print("Reynolds (Re_in): \t\t\t%.2E"%(Re_in))
print("Resistência de superfície (Lambda): \t%.4f"%(Lambda))
print("Parâmetro Geométrico real (A_eq): \t%.2f"%(A_eq))
print("Coef. de área livre real (Phi_eq): \t%.4f"%(Phi_eq))
print("Espessura da camada líquida (t_as_w): \t%.4f \tmm"%(t_as_w))
print("Perda de momento angular (K): \t\t%.4f"%(K))
print("Redução de área do estágio (C): \t%.4f"%(C))
print("Coef. de perda (Ksi_in): \t\t%.4f"%(Ksi_in))
print("Perda nos canais (Ksi): \t\t%.4f"%(Ksi))
print("Coeficiente de descarga real (Mu_w): \t%.4f"%(Mu_w))
print("Coeficinete de descarga total (Mu_i): \t%.4f"%(Mu_i))
print("Vazão mássica real (m_r): \t\t%.1f \tg/s"%(m_r))

##### Cálulo do ângulo de spray real #####

Alpha2 = 2*180/np.pi*np.arcsin(2*Mu_i*A_eq/((1+np.sqrt(1-Phi_i))*np.sqrt(1-Ksi*Mu_i**2*A_i**2/C**2)))


print("Ângulo de Abertura do Spray (Alpha2): \t%.1f \tdeg"%(Alpha2))

##### Cálculo da não-uniformidade #####

I = 23.7/((R_in/R_n)**2.7*N_in**1.34*Phi_eq**1.1*(L_cv/D_cv)**0.15)

print("\n\nNão Uniformidade Esperada = \t%.2f %%"%(I))


##### Raios do fluido #####

r_mn = R_n*np.sqrt(1-Phi_eq)
    
t_fluid = R_n-r_mn
    
r_mk = r_mn * np.sqrt(2*(1-Phi_eq)/(2-Phi_eq))
print("r_mk = ",r_mk)



print("\n\n\nCurva de vazão:\n")


##### plot da curva de vazão ####

Mass_ideal = []
Mass_real = []
dPressure = np.arange(0.1,0.95,0.01)
Alfa = []

for Delta_P in dPressure:
    ### Ideal Injector Parameters  ###
    
    A_i = calc_A (R_in, R_n, N_in, r_in, 0)
    
    ## Phi_i calculation ##
     # Start #
    
    Phi_i = calc_Phi(A_i)
    
     # End #
     
    #Coeficiente de Descarga ideal
    Mu_i = Phi_i*np.sqrt(Phi_i/(2-Phi_i))
    
    #Espessura de camada líquida ideal
    t_as_i = R_n*(1-np.sqrt(1-Phi_i)) #mm
    
    if(t_as_i <= 0):
        print("ERROR: invalid t_as value")
        
    elif(t_as_i <= 0.1):
        print("Caution! Low value of t_as")
        
    #Vazão mássica ideal
    m_i = Mu_i*A_n*1e-6*np.sqrt(2*Delta_P*1e6*Rho_p)*1e3 #vazão mássica ideal [g/s]
    
    #####   Cálculo do Estágio Real   #####
    
    
    #Velocidade Tangencial
    W = m_i*1e-3/(N_in*np.pi*(r_in*1e-3)**2*Rho_p) #Tangential Velocity
    
    
    
    #Número de Reynolds nos canais tangenciais
    Re_in = W*d_in*1e-3/(Din_visc_p/Rho_p)
    
    
    
    if(Re_in < 2600): print("ERROR: Low Reynolds number (%.2E < 2600) [\u0394P = %.2f MPa]"%(Re_in,Delta_P))
    elif(Re_in > 5e4): print("Caution: Reynolds number over 5e4")
    
    #Coeficiente de resistência de superfície
    Lambda = 0.3164*Re_in**(-0.25)
    
    
    
    #Parâmetro Geométrico equivalente
    A_eq = calc_A (R_in, R_n, N_in, r_in, Lambda) 
    
    
    
    ##### Cálculo do Coeficiente de área livre real
    

    Phi_eq = calc_Phi(A_eq)
    
    #####
    
    #Coeficiente de Descarga
    
    Mu_w = Phi_eq*np.sqrt(Phi_eq)/(np.sqrt(2-Phi_eq))
    
    
    #Espessura de camada líquida real
    t_as_w = R_n*(1-np.sqrt(1-Phi_eq)) #mm
    
    
    
    if(t_as_w <= 0):
        print("ERROR: invalid fluid layer thickness value")
        
    elif(t_as_w <= 0.1):
        print("Caution! Low fluid layer thickness")
    
    
    #Coeficiente de perda de momento angular
    K = A_eq/A_i
    #if(K>1): print("ERROR: K>1")
        
    
    
    #Coeficiente de redução de área do estágio
    C = R_in/R_n
    
    
    ##### Cálculo da perda nos canais tangenciais #####
    
    
    Alpha_in = 90 - 180/np.pi*np.arctan(R_cv/L_in)
    
    #Coeficiente de perda
    Ksi_in = -0.0204*(Fillet/d_in)**4-0.026*(Fillet/d_in)**3+0.3587*(Fillet/d_in)**2-0.6995*(Fillet/d_in)+0.5002
    #Ksi_in = -1*Alpha_in/150 + 1.1
    
    
    #Fator de perda nos canais tangenciais
    Ksi = Ksi_in + Lambda*L_in/(2*r_in)
    
     #coeficiente de descarga total
    Mu_i = 1/np.sqrt(1/Phi_eq**2+A_i**2/(1-Phi_eq)+(Ksi*N_in)*A_i**2/C**2)
    #Mu_i = Mu_w/(np.sqrt(1+Mu_w**2*Ksi*A_i**2/C**2))
    
    #Vazão mássica real
    m_r = Mu_i*np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P0*1e6)*1000 #np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P*1e6)/np.sqrt((2-Phi_eq)/Phi_eq**3+Ksi*A_i**2/R_in**2)*1e3
    #m_r = Mu_i*np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P*1e6)*1e3
    
    
    Mass_ideal.append(m_i)
    
    Mass_real.append(m_r)
    
    Alfa_seq = 2*180/np.pi*np.arcsin(2*Mu_i*A_eq/((1+np.sqrt(1-Phi_eq))*np.sqrt(1-Ksi*Mu_i**2*A_i**2/C**2)))#2*180/np.pi*np.arcsin(2*Mu_i*A_i*K/(1+np.sqrt(1+Phi_eq)*np.sqrt(1-Ksi*Mu_i**2*A_i**2)/C**2)) if -1 <= 2*Mu_i*A_i*K/(1+np.sqrt(1+Phi_eq)*np.sqrt(1-Ksi*Mu_i**2*A_i**2)/C**2) <= 1 else 0
    
    Alfa.append(Alfa_seq)
    
M_target = []

for i in dPressure:
    M_target.append(m_target)

for i in range(len(dPressure)):
    if(abs(dPressure[i] - Delta_P0) < 0.005):
        print("\nm_bar_i = %.2f g/s"%(Mass_ideal[i]))
        print("\nm_bar = %.2f g/s"%(Mass_real[i]))
        print("%%target = %.2f %%"%(Mass_real[i]/m_target*100))



plt.figure(1)
line01 = plt.plot(dPressure, Mass_ideal, 'k--', markersize=1, label = "ṁ fluido ideal")
line02 = plt.plot(dPressure, Mass_real, 'k-', markersize=1, label = "ṁ fluido viscoso")
line03 = plt.plot(dPressure, M_target, 'k-.', markersize=1, label = "ṁ alvo")
plt.legend()
plt.ylabel("ṁ [g/s]")
plt.xlabel("ΔP [MPa]")
plt.xticks(np.arange(0.1,1.0+0.1,0.1))
plt.grid()


'''
plt.figure(2)
plt.plot(dPressure, Alfa, markersize=1)
plt.xticks(np.arange(0.1,1.0+0.1,0.1))
plt.grid()
'''




##### ##### Cálculos de Tamanho de Gota ##### #####

#Radcliffe [16]
#


#Jasuja [15]
#


#Babu et al. [37]
#


#Jones [10]


#Lefebvre [13]


#Wang and Lefebvre [14]




