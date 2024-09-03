# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 07:41:29 2020

@author: A. Goulart
"""

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt

from M2_FileReader import FileReader_2


class Method_2:
    def __init__(self, foldername, filename):
        self.it_lim = 1000
        
        self.fr = FileReader_2()
        self.fr.setFolderName(foldername)
        self.fr.setFileName(filename)
        
        ###Relações úteis###
        self.rad_to_deg = 180/np.pi
        self.deg_to_rad = np.pi/180
    
    
    def calc_Phi(self,A):
        self.Phi_0 = 0.1 #initial value
        
        def f_Phi(x):
            return (1-x)*np.sqrt(2)/x**1.5-A
        
        self.Phi = opt.fsolve(f_Phi, self.Phi_0)
        
        if(len(self.Phi) == 1):
            self.Phi = self.Phi[0]
            return self.Phi
            
        elif(len(self.Phi) == 0):
            print("ERROR: no Phi value found")
            
        else:
            print("ERROR: Phi multiple values")
    
    
    
    def calc_A (self,R_in, R_n, N_in, r_in, Lambda):
        return R_in*1e-3*R_n*1e-3/(N_in*(r_in*1e-3)**2+Lambda/2*R_in*1e-3*(R_in*1e-3-R_n*1e-3))
    
    
    def run_M2(self):
        
        
        
        self.fr.read()
        
        
        #SIG - Swirl Injector from Geometry
        print("\n Method M2 \n")

        ##### ##### ##### ##### ##### ##### #####
        ##### #####  ST1 Calculation  ##### #####
        ##### ##### ##### ##### ##### ##### #####
        
        
        print("m_target_1: \t%.2f \tg/s\n"%(np.round(self.fr.m_1,2)))
        


        print("\n ST1 - Geometria do injetor:")
 
        print("Número de canais de entrada (N_1): \t\t\t\t%.1f"%(self.fr.N_1))

        print("Diâmetro dos canais de entrada (d_in_1): \t\t%.2f \tmm"%(self.fr.d_in_1))

        print("Comprimento dos canais de entrada (L_in_1): \t%.2f \tmm"%(self.fr.L_in_1))

        print("Diâmetro da Câmara de Vórtice (d_s_1): \t\t\t%.2f \tmm"%(self.fr.d_s_1))

        print("Comprimento da Câmara de Vórtice (L_s_1): \t\t%.2f \tmm"%(self.fr.L_s_1))
        
        print("Diâmetro do Orifício de saída (D_n): \t\t\t%.2f \tmm"%(self.fr.d_n_1))
        
        print("Raio de Entrada dos canais (R_0_1): \t\t\t%.2f \tmm"%(self.fr.R_0_1))
        
        print("\nDelta_P0_1: \t\t%.2f bar"%(self.fr.Delta_P0_1))
        
        
        ##### ##### ##### ##### ##### ###### #####
        ##### ST1 - Cálculo do Estágio Ideal #####
        ##### ##### ##### ##### ##### ###### #####
        
        print("\n\n Dados do Estágio ideal:")
        
        ### Ideal Injector Parameters  ###
        
        self.A_i = self.calc_A(self.fr.R_0_1, self.fr.d_n_1/2, self.fr.N_1, self.fr.d_in_1/2, 0) 
         
        
        self.Phi_i = self.calc_Phi(self.A_i)
        
         
         
        #Coeficiente de Descarga ideal
        self.Mu_i = self.Phi_i*np.sqrt(self.Phi_i/(2-self.Phi_i))
        
        #Espessura de camada líquida ideal
        self.t_as_i = self.fr.R_0_1*(1-np.sqrt(1-self.Phi_i)) #mm
        
        if(self.t_as_i <= 0):
            print("ERROR: invalid t_as value")
            
        elif(self.t_as_i <= 0.1):
            print("Caution! Low value of t_as")
            
        #Vazão mássica ideal
        self.m_i = self.Mu_i*(np.pi/4*self.fr.d_n_1**2)*1e-6*np.sqrt(2*self.fr.Delta_P0_1*1e5*self.fr.Rho_1)*1e3 #vazão mássica ideal [g/s]
        
        
        #Coeficiente de redução de área do estágio
        self.C = self.fr.R_0_1/(self.fr.d_n_1/2)
        
        
        print("Parâmetro Geométrico ideal (A_i): \t\t%.2f"%(self.A_i))
        print("Coef. de área livre ideal (Phi_i): \t\t%.4f"%(self.Phi_i))
        print("Coeficiente de descarga ideal(Mu_i): \t%.4f"%(self.Mu_i))
        print("Espessura da camada líquida (t_as_i): \t%.4f \tmm"%(self.t_as_i))
        print("Vazão mássica ideal (m_i): \t\t%.1f \tg/s"%(self.m_i))
        
        ##### Cálulo do ângulo de spray ideal #####
        
        self.Alpha2_i = 2*180/np.pi*np.arcsin(2*self.Mu_i*self.A_i/((1+np.sqrt(1-self.Phi_i))*np.sqrt(1-0*self.Mu_i**2*self.A_i**2/self.C**2)))
        
        print("Ângulo de Abertura do Spray (Alpha2_i): \t%.1f \tdeg"%(self.Alpha2_i))
        
        ##### ##### ##### ##### ##### ##### #####
        ##### ST1 - Cálculo do Estágio Real #####
        ##### ##### ##### ##### ##### ##### #####
        
        print("\n Dados do Estágio real:")
        
        
        #Velocidade Tangencial
        self.W = self.m_i*1e-3/(self.fr.N_1*(np.pi/4*self.fr.d_in_1**2)*1e-6*self.fr.Rho_1) #Tangential Velocity
        
        
        
        #Número de Reynolds nos canais tangenciais
        self.Re_in = self.fr.Rho_1*self.W*self.fr.d_in_1*1e-3/self.fr.Din_visc_1
        
        
        
        if(self.Re_in < 2600): print("ERROR: Low Reynolds number (%.2E < 2600)"%(self.Re_in))
        elif(self.Re_in > 5e4): print("Caution: Reynolds number over 5e4")
        
        #Coeficiente de resistência de superfície
        self.Lambda = 0.3164*self.Re_in**(-0.25) #Bazarov
        #Lambda = 10**(25.8/(np.log10(Re_in))**2.58-2) #Bayvel
        
        
        #Parâmetro Geométrico equivalente
        self.A_w = self.calc_A(self.fr.R_0_1, self.fr.d_n_1/2, self.fr.N_1, self.fr.d_in_1/2, self.Lambda) 
        
        
        
        ##### Cálculo do Coeficiente de área livre real
        
        self.Phi_w = self.calc_Phi(self.A_w)
        
        
        #####
        
        #Coeficiente de Descarga real
            
        self.Mu_w = self.Phi_w**1.5/(np.sqrt(2-self.Phi_w))
        
        
        #Espessura de camada líquida real
        self.t_as_w = self.fr.d_n_1/2*(1-np.sqrt(1-self.Phi_w)) #mm
        
        
        
        if(self.t_as_w <= 0):
            print("ERROR: invalid fluid layer thickness value")
            
        elif(self.t_as_w <= 0.1):
            print("Caution! Low fluid layer thickness")
        
        
        #Coeficiente de perda de momento angular
        self.K = self.A_w/self.A_i
        
        
        
        ##### Cálculo da perda nos canais tangenciais #####
        
        self.Alpha_in = 90 - 180/np.pi*np.arctan(self.fr.d_s_1/2/self.fr.L_in_1)
        
        #Coeficiente de perda
        self.Ksi_in = -1*self.Alpha_in/150 + 1.1   #-0.0204*(d_bx/d_in)**4-0.026*(d_bx/d_in)**3+0.3587*(d_bx/d_in)**2-0.6995*(d_bx/d_in)+0.5002
        
        
        
        #Fator de perda nos canais tangenciais
        self.Ksi = self.Ksi_in + self.Lambda*self.fr.L_in_1/(self.fr.d_in_1)
        
        
        
        #Vazão mássica real
        self.m_r = np.pi*(self.fr.d_n_1/2*1e-3)**2*np.sqrt(2*self.fr.Rho_1*self.fr.Delta_P0_1*1e5)/np.sqrt((2-self.Phi_w)/self.Phi_w**3+self.Ksi*self.A_i**2/(self.fr.R_0_1)**2)*1e3   #np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P0*1e5)/np.sqrt(1/Phi_w**2+A_w**2/(1-Phi_w)+(Ksi*N_in)*A_i**2/C**2)*1000
        
        #coeficiente de descarga total
        self.Mu_t = self.Mu_w/(np.sqrt(1+self.Mu_w**2*self.Ksi*self.A_i**2/self.C**2))
        
        print("Velocidade tangencial (W): \t\t\t\t%.2f \tm/s"%(self.W))
        print("Reynolds (Re_in): \t\t\t\t\t\t%.2E"%(self.Re_in))
        print("Resistência de superfície (Lambda): \t%.4f"%(self.Lambda))
        print("Parâmetro Geométrico real (A_w): \t\t%.2f"%(self.A_w))
        print("Coef. de área livre real (Phi_w): \t\t%.4f"%(self.Phi_w))
        print("Espessura da camada líquida (t_as_w): \t%.4f \tmm"%(self.t_as_w))
        print("Perda de momento angular (K): \t\t\t%.4f"%(self.K))
        print("Redução de área do estágio (C): \t\t%.4f"%(self.C))
        print("Coef. de perda (Ksi_in): \t\t\t\t%.4f"%(self.Ksi_in))
        print("Perda nos canais (Ksi): \t\t\t\t%.4f"%(self.Ksi))
        print("Coeficiente de descarga real (Mu_w): \t%.4f"%(self.Mu_w))
        print("Coeficinete de descarga total (Mu_t): \t%.4f"%(self.Mu_t))
        print("Vazão mássica real (m_r): \t\t\t\t%.1f \tg/s"%(self.m_r))
        
        
        print("\nDiâmetro da saída (R_out_2): \t\t\t%.2f mm"%(self.fr.R_out_2))
        
        ##### Cálulo do ângulo de spray real #####
        
        self.Alpha2 = 2*180/np.pi*np.arcsin(2*self.Mu_i*self.A_w/((1+np.sqrt(1-self.Phi_i))*np.sqrt(1-self.Ksi*self.Mu_i**2*self.A_i**2/self.C**2)))
        
        print("Ângulo de Abertura do Spray (Alpha2): \t%.1f \tdeg"%(self.Alpha2))
        
        
        ##### Cálculo da não-uniformidade #####
        
        self.I = 23.7/((self.fr.R_0_1/(self.fr.d_n_1/2))**2.7*self.fr.N_1**1.34*self.Phi_w**1.1*(self.fr.L_s_1/self.fr.d_s_1)**0.15)
        
        print("\n\nNão Uniformidade Esperada = \t%.2f %%"%(self.I))
        
        
        
        
        ##### ##### ##### ##### ##### ##### ##### ##### #####
        ##### ##### ST1 - plot da curva de vazão  ##### #####
        ##### ##### ##### ##### ##### ##### ##### ##### #####
        
        print("\n\n\nCurva de vazão:\n")
        
        self.Mass_ideal = []
        self.Alpha2_ideal = []
        self.Mass_real = []
        self.Alpha2_real = []
        self.Re_plot = []
        self.Ksi_plot = []
        self.dPressure = np.arange(0.50,9.90,0.01) #bar

        
        for self.Delta_P in self.dPressure:
            ### Ideal Injector Parameters  ###
            
            self.A_i = self.calc_A (self.fr.R_0_1, self.fr.d_n_1/2, self.fr.N_1, (self.fr.d_in_1/2), 0)
            
            ## Phi_i calculation ##
                        
            self.Phi_i = self.calc_Phi(self.A_i)
            
             
            #Coeficiente de Descarga ideal
            self.Mu_i = self.Phi_i*np.sqrt(self.Phi_i/(2-self.Phi_i))
            
            #Espessura de camada líquida ideal
            self.t_as_i = self.fr.d_n_1/2*(1-np.sqrt(1-self.Phi_i)) #mm
            
            if(self.t_as_i <= 0):
                print("ERROR: invalid t_as value")
                
            elif(self.t_as_i <= 0.1):
                print("Caution! Low value of t_as")
                
            #Vazão mássica ideal
            self.m_i = self.Mu_i*(np.pi/4*self.fr.d_n_1**2)*1e-6*np.sqrt(2*self.Delta_P*1e5*self.fr.Rho_1)*1e3 #vazão mássica ideal [g/s]
            
            #####   Cálculo do Estágio Real   #####
                        
            #Velocidade Tangencial
            self.W = self.m_i*1e-3/(self.fr.N_1*(np.pi/4*self.fr.d_in_1**2)*1e-6*self.fr.Rho_1) #Tangential Velocity

            
            
            #Número de Reynolds nos canais tangenciais
            self.Re_in = self.fr.Rho_1*self.W*self.fr.d_in_1*1e-3/self.fr.Din_visc_1
            
            

            if(self.Re_in < 2600): print("ERROR: Low Reynolds number (%.2E < 2600) [\u0394P = %.2f bar]"%(self.Re_in,self.Delta_P))
            elif(self.Re_in > 5e4): print("Caution: Reynolds number over 5e4")

            #Coeficiente de resistência de superfície
            self.Lambda = 0.3164*self.Re_in**(-0.25) #Bazarov
            #Lambda = 10**(25.8/(np.log10(Re_in))**2.58-2) #Bayvel
            
            
            
            #Parâmetro Geométrico equivalente
            self.A_w = self.calc_A (self.fr.R_0_1, (self.fr.d_n_1/2), self.fr.N_1, (self.fr.d_in_1/2), self.Lambda) 
            
            
            
            ##### Cálculo do Coeficiente de área livre real
            
        
            self.Phi_w = self.calc_Phi(self.A_w)
            
            #####
            
            #Coeficiente de Descarga
            
            self.Mu_w = self.Phi_w*np.sqrt(self.Phi_w)/(np.sqrt(2-self.Phi_w))
            
            
            #Espessura de camada líquida real
            self.t_as_w = self.fr.d_n_1/2*(1-np.sqrt(1-self.Phi_w)) #mm
            
            
            
            if(self.t_as_w <= 0):
                print("ERROR: invalid fluid layer thickness value")
                
            elif(self.t_as_w <= 0.1):
                print("Caution! Low fluid layer thickness")
            
            
            #Coeficiente de perda de momento angular
            self.K = self.A_w/self.A_i
            
            
            
            #Coeficiente de redução de área do estágio
            self.C = self.fr.R_0_1/(self.fr.d_n_1/2)
            
            
            ##### Cálculo da perda nos canais tangenciais #####
            
            
            self.Alpha_in = 90 - 180/np.pi*np.arctan(self.fr.d_s_1/2/self.fr.L_in_1)
            
            #Coeficiente de perda
            self.Ksi_in = -1*self.Alpha_in/150 + 1.1#-0.0204*(d_bx/d_in)**4-0.026*(d_bx/d_in)**3+0.3587*(d_bx/d_in)**2-0.6995*(d_bx/d_in)+0.5002
            
            
            
            #Fator de perda nos canais tangenciais
            self.Ksi = self.Ksi_in + self.Lambda*self.fr.L_in_1/(self.fr.d_in_1)
            
            
            
            #Vazão mássica real
            self.m_r = np.pi/4*(self.fr.d_n_1*1e-3)**2*np.sqrt(2*self.fr.Rho_1*self.Delta_P*1e5)/np.sqrt((2-self.Phi_w)/self.Phi_w**3+self.Ksi*self.A_i**2/self.fr.R_0_1**2)*1e3  #np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P*1e5)/np.sqrt(1/Phi_w**2+A_i**2/(1-Phi_w)+(Ksi*N_in)*A_i**2/C**2)*1000
            
            #coeficiente de descarga total
            self.Mu_t = self.Mu_w/(np.sqrt(1+self.Mu_w**2*self.Ksi*self.A_i**2/self.C**2))
            
            
            self.Mass_ideal.append(self.m_i)
        
            
            
            self.Mass_real.append(self.m_r)
            
            self.Alpha2_i = 2*180/np.pi*np.arcsin(2*self.Mu_i*self.A_i/((1+np.sqrt(1-self.Phi_i))*np.sqrt(1-0*self.Mu_i**2*self.A_i**2/self.C**2)))
            
            self.Alpha2_w = 2*180/np.pi*np.arcsin(2*self.Mu_i*self.A_w/((1+np.sqrt(1-self.Phi_i))*np.sqrt(1-self.Ksi*self.Mu_i**2*self.A_i**2/self.C**2)))
            
            #self.Alpha2_w = 2*180/np.pi*np.arcsin(2*self.Mu_t*self.A_i*self.K/(1+np.sqrt(1+self.Phi_w)*np.sqrt(1-self.Ksi*self.Mu_t**2*self.A_i**2)/self.C**2)) if -1 <= 2*self.Mu_t*self.A_i*self.K/(1+np.sqrt(1+self.Phi_w)*np.sqrt(1-self.Ksi*self.Mu_t**2*self.A_i**2)/self.C**2) <= 1 else 0
            
            self.Alpha2_ideal.append(self.Alpha2_i)
            
            self.Alpha2_real.append(self.Alpha2_w)
            
            self.Re_plot.append(self.Re_in)
            
            self.Ksi_plot.append(self.Ksi)
            

            
            self.M_target = []
            
            for self.i in self.dPressure:
                self.M_target.append(self.fr.m_1)
            
        for self.i in range(len(self.dPressure)):
            if(abs(self.dPressure[self.i] - self.fr.Delta_P0_1) < 0.005):

                print("\nm_bar_i = %.2f g/s"%(self.Mass_ideal[self.i]))
                print("\nm_bar = %.2f g/s"%(self.Mass_real[self.i]))
                print("%%target = %.2f %%"%(self.Mass_real[self.i]/self.fr.m_1*100))
            
            
        fig, ax11 = plt.subplots(figsize=(8, 6))
        ax11.grid()
        ax11.set_title("ST1 - ṁ x ΔP") 
        ax11.plot(self.dPressure, self.Mass_ideal, 'k--', markersize=1, label = "ṁ fluido ideal", color = '0.5')
        ax11.plot(self.dPressure, self.Mass_real, 'k:', markersize=3, label = "ṁ fluido viscoso")
        ax11.plot(self.dPressure, self.M_target, 'k-.', markersize=1, label = "ṁ alvo")
        ax11.legend()
        ax11.set_ylabel("ṁ [g/s]")
        ax11.set_xlabel("ΔP [bar]")
        ax11.set_xticks(np.arange(0,10.0+0.5,0.5))
        ax11.set_xlim(0.0,10.0)
        #ax11.grid()
            
        fig, ax12 = plt.subplots(figsize=(8, 6))
        ax12.grid()
        ax12.set_title("ST1 - Re x ΔP") 
        ax12.plot(self.dPressure, self.Re_plot, 'k-', markersize=1)
        ax12.set_ylabel("Re [ad.]")
        ax12.set_xlabel("ΔP [bar]")
        ax12.set_xticks(np.arange(0,10.0+0.5,0.5))
        ax12.set_xlim(0.0,10.0)
        
        fig, ax13 = plt.subplots(figsize=(8, 6))
        ax13.grid()
        ax13.set_title("ST1 - Ksi x ΔP") 
        ax13.plot(self.dPressure, self.Ksi_plot, 'k-', markersize=1)
        ax13.set_ylabel("Ksi [ad.]")
        ax13.set_xlabel("ΔP [bar]")
        ax13.set_xticks(np.arange(0,10.0+0.5,0.5))
        ax13.set_xlim(0.0,10.0)
        
        fig, ax14 = plt.subplots(figsize=(8, 6))
        ax14.grid()
        ax14.set_title("ST1 - 2α x ΔP") 
        ax14.plot(self.dPressure, self.Alpha2_ideal, 'k--', markersize=1, label = "2α fluido ideal", color = '0.5')
        ax14.plot(self.dPressure, self.Alpha2_real, 'k:', markersize=1, label = "2α fluido viscoso")
        ax14.legend()
        ax14.set_ylabel("2α [deg]")
        ax14.set_xlabel("ΔP [bar]")
        ax14.set_xticks(np.arange(0,10.0+0.5,0.5))
        ax14.set_xlim(0.0,10.0)
        
        ##### ##### ##### ##### ##### ##### #####
        ##### #####  ST2 Calculation  ##### #####
        ##### ##### ##### ##### ##### ##### #####
            
        
        print("m_target_2: \t%.2f \tg/s\n"%(np.round(self.fr.m_2,2)))
        


        print("\n ST2 - Geometria do injetor:")
 
        print("Número de canais de entrada (N_2): \t\t\t\t%.1f"%(self.fr.N_2))

        print("Diâmetro dos canais de entrada (d_in_2): \t\t%.2f \tmm"%(self.fr.d_in_2))

        print("Comprimento dos canais de entrada (L_in_2): \t%.2f \tmm"%(self.fr.L_in_2))

        print("Diâmetro da Câmara de Vórtice (d_s_2): \t\t\t%.2f \tmm"%(self.fr.d_s_2))

        print("Comprimento da Câmara de Vórtice (L_s_2): \t\t%.2f \tmm"%(self.fr.L_s_2))
        
        print("Diâmetro do Orifício de saída (D_n): \t\t\t%.2f \tmm"%(self.fr.d_n_2))
        
        print("Raio de Entrada dos canais (R_0_2): \t\t\t%.2f \tmm"%(self.fr.R_0_2))
        
        
        print("\nDelta_P0_2: \t\t%.2f bar"%(self.fr.Delta_P0_2))
        
        
        ##### ##### ##### ##### ##### ###### #####
        ##### ST2 - Cálculo do Estágio Ideal #####
        ##### ##### ##### ##### ##### ###### #####
        
        print("\n\n Dados do Estágio ideal:")
        
        ### Ideal Injector Parameters  ###
        
        self.A_i = self.calc_A(self.fr.R_0_2, self.fr.d_n_2/2, self.fr.N_2, self.fr.d_in_2/2, 0) 
         
        
        self.Phi_i = self.calc_Phi(self.A_i)
        
         
         
        #Coeficiente de Descarga ideal
        self.Mu_i = self.Phi_i*np.sqrt(self.Phi_i/(2-self.Phi_i))
        
        #Espessura de camada líquida ideal
        self.t_as_i = self.fr.R_0_2*(1-np.sqrt(1-self.Phi_i)) #mm
        
        if(self.t_as_i <= 0):
            print("ERROR: invalid t_as value")
            
        elif(self.t_as_i <= 0.1):
            print("Caution! Low value of t_as")
            
        #Vazão mássica ideal
        self.m_i = self.Mu_i*(np.pi/4*self.fr.d_n_2**2)*1e-6*np.sqrt(2*self.fr.Delta_P0_2*1e5*self.fr.Rho_2)*1e3 #vazão mássica ideal [g/s]
        
        
        #Coeficiente de redução de área do estágio
        self.C = self.fr.R_0_2/(self.fr.d_n_2/2)
        
        
        print("Parâmetro Geométrico ideal (A_i): \t\t%.2f"%(self.A_i))
        print("Coef. de área livre ideal (Phi_i): \t\t%.4f"%(self.Phi_i))
        print("Coeficiente de descarga ideal(Mu_i): \t%.4f"%(self.Mu_i))
        print("Espessura da camada líquida (t_as_i): \t%.4f \tmm"%(self.t_as_i))
        print("Vazão mássica ideal (m_i): \t\t%.1f \tg/s"%(self.m_i))
        
        ##### Cálulo do ângulo de spray ideal #####
        
        self.Alpha2 = 2*180/np.pi*np.arcsin(2*self.Mu_i*self.A_i/((1+np.sqrt(1-self.Phi_i))*np.sqrt(1-0*self.Mu_i**2*self.A_i**2/self.C**2)))
        
        print("Ângulo de Abertura do Spray (Alpha2): \t%.1f \tdeg"%(self.Alpha2))
        
        ##### ##### ##### ##### ##### ##### #####
        ##### ST2 - Cálculo do Estágio Real #####
        ##### ##### ##### ##### ##### ##### #####
        
        print("\n Dados do Estágio real:")
        
        
        #Velocidade Tangencial
        self.W = self.m_i*1e-3/(self.fr.N_2*(np.pi/4*self.fr.d_in_2**2)*1e-6*self.fr.Rho_2) #Tangential Velocity
        
        
        
        #Número de Reynolds nos canais tangenciais
        self.Re_in = self.fr.Rho_2*self.W*self.fr.d_in_2*1e-3/self.fr.Din_visc_2
        
        
        
        if(self.Re_in < 2600): print("ERROR: Low Reynolds number (%.2E < 2600)"%(self.Re_in))
        elif(self.Re_in > 5e4): print("Caution: Reynolds number over 5e4")
        
        #Coeficiente de resistência de superfície
        self.Lambda = 0.3164*self.Re_in**(-0.25) #Bazarov
        #Lambda = 10**(25.8/(np.log10(Re_in))**2.58-2) #Bayvel
        
        
        #Parâmetro Geométrico equivalente
        self.A_w = self.calc_A(self.fr.R_0_2, (self.fr.d_n_2/2), self.fr.N_2, self.fr.d_in_2/2, self.Lambda) 
        
        
        
        ##### Cálculo do Coeficiente de área livre real
        
        self.Phi_w = self.calc_Phi(self.A_w)
        
        
        #####
        
        #Coeficiente de Descarga
            
        self.Mu_w = self.Phi_w**1.5/(np.sqrt(2-self.Phi_w))
        
        
        #Espessura de camada líquida real
        self.t_as_w = self.fr.d_n_2/2*(1-np.sqrt(1-self.Phi_w)) #mm
        
        
        
        if(self.t_as_w <= 0):
            print("ERROR: invalid fluid layer thickness value")
            
        elif(self.t_as_w <= 0.1):
            print("Caution! Low fluid layer thickness")
        
        
        #Coeficiente de perda de momento angular
        self.K = self.A_w/self.A_i
        
        
        
        ##### Cálculo da perda nos canais tangenciais #####
        
        self.Alpha_in = 90 - 180/np.pi*np.arctan(self.fr.d_s_2/2/self.fr.L_in_2)
        
        #Coeficiente de perda
        self.Ksi_in = -1*self.Alpha_in/150 + 1.1   #-0.0204*(d_bx/d_in)**4-0.026*(d_bx/d_in)**3+0.3587*(d_bx/d_in)**2-0.6995*(d_bx/d_in)+0.5002
        
        
        
        #Fator de perda nos canais tangenciais
        self.Ksi = self.Ksi_in + self.Lambda*self.fr.L_in_2/(self.fr.d_in_2)
        
        
        
        #Vazão mássica real
        self.m_r = np.pi/4*(self.fr.d_n_2*1e-3)**2*np.sqrt(2*self.fr.Rho_2*self.fr.Delta_P0_2*1e5)/np.sqrt((2-self.Phi_w)/self.Phi_w**3+self.Ksi*self.A_i**2/(self.fr.R_0_2)**2)*1e3   #np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P0*1e5)/np.sqrt(1/Phi_w**2+A_w**2/(1-Phi_w)+(Ksi*N_in)*A_i**2/C**2)*1000
        
        #coeficiente de descarga total
        self.Mu_t = self.Mu_w/(np.sqrt(1+self.Mu_w**2*self.Ksi*self.A_i**2/self.C**2))
        
        print("Velocidade tangencial (W): \t\t\t\t%.2f \tm/s"%(self.W))
        print("Reynolds (Re_in): \t\t\t\t\t\t%.2E"%(self.Re_in))
        print("Resistência de superfície (Lambda): \t%.4f"%(self.Lambda))
        print("Parâmetro Geométrico real (A_w): \t\t%.2f"%(self.A_w))
        print("Coef. de área livre real (Phi_w): \t\t%.4f"%(self.Phi_w))
        print("Espessura da camada líquida (t_as_w): \t%.4f \tmm"%(self.t_as_w))
        print("Perda de momento angular (K): \t\t\t%.4f"%(self.K))
        print("Redução de área do estágio (C): \t\t%.4f"%(self.C))
        print("Coef. de perda (Ksi_in): \t\t\t\t%.4f"%(self.Ksi_in))
        print("Perda nos canais (Ksi): \t\t\t\t%.4f"%(self.Ksi))
        print("Coeficiente de descarga real (Mu_w): \t%.4f"%(self.Mu_w))
        print("Coeficinete de descarga total (Mu_t): \t%.4f"%(self.Mu_t))
        print("Vazão mássica real (m_r): \t\t\t\t%.1f \tg/s"%(self.m_r))
        #print("Test1 Mu_w: %.2f"%(np.pi/4*(self.fr.d_n_2*1e-3)**2*self.Mu_w*np.sqrt(2*self.fr.Rho_2*self.fr.Delta_P0_2*1e5)*1e3))
        #print("Test1 Mu_t: %.2f"%(np.pi/4*(self.fr.d_n_2*1e-3)**2*self.Mu_t*np.sqrt(2*self.fr.Rho_2*self.fr.Delta_P0_2*1e5)*1e3))
        
        ##### Cálulo do ângulo de spray real #####
        
        self.Alpha2 = 2*180/np.pi*np.arcsin(2*self.Mu_i*self.A_w/((1+np.sqrt(1-self.Phi_i))*np.sqrt(1-self.Ksi*self.Mu_i**2*self.A_i**2/self.C**2)))
        

        
        #print("Test1: ", 2*self.Mu_i*self.A_i/((1+np.sqrt(1-self.Phi_i))*np.sqrt(1-0*self.Mu_i**2*self.A_i**2/self.C**2)))
        #print("Test2: ", 2*self.Mu_i*self.A_w/((1+np.sqrt(1-self.Phi_i))*np.sqrt(1-self.Ksi*self.Mu_i**2*self.A_i**2/self.C**2)))
        
        print("Ângulo de Abertura do Spray (Alpha2): \t%.1f \tdeg"%(self.Alpha2))
        
        
        self.a = 2*(1-self.Phi_w)**2/(2-self.Phi_w)
        
        self.Alpha2_out = np.arctan(self.a/((self.fr.R_out_2/(self.fr.d_n_2/2))**2-self.a))
        
        print("Ângulo de Abertura do Spray modificado (Alpha2_out): \t%.1f \tdeg"%(self.Alpha2_out))
        
        ##### Cálculo da não-uniformidade #####
        
        self.I = 23.7/((self.fr.R_0_2/(self.fr.d_n_2/2))**2.7*self.fr.N_2**1.34*self.Phi_w**1.1*(self.fr.L_s_2/self.fr.d_s_2)**0.15)
        
        print("\n\nNão Uniformidade Esperada = \t%.2f %%"%(self.I))
        
        
        
        
        ##### ##### ##### ##### ##### ##### ##### ##### #####
        ##### ##### ST2 - plot da curva de vazão  ##### #####
        ##### ##### ##### ##### ##### ##### ##### ##### #####
        
        print("\n\n\nCurva de vazão:\n")
        
        self.Mass_ideal = []
        self.Alpha2_ideal = []
        self.Mass_real = []
        self.Alpha2_real = []
        self.Re_plot = []
        self.Ksi_plot = []
        self.test = []
        self.dPressure = np.arange(0.50,9.90,0.01) #bar

        
        for self.Delta_P in self.dPressure:
            ### Ideal Injector Parameters  ###
            
            self.A_i = self.calc_A (self.fr.R_0_2, self.fr.d_n_2/2, self.fr.N_2, (self.fr.d_in_2/2), 0)
            
            ## Phi_i calculation ##
                        
            self.Phi_i = self.calc_Phi(self.A_i)
            
             
            #Coeficiente de Descarga ideal
            self.Mu_i = self.Phi_i*np.sqrt(self.Phi_i/(2-self.Phi_i))
            
            #Espessura de camada líquida ideal
            self.t_as_i = self.fr.d_n_2/2*(1-np.sqrt(1-self.Phi_i)) #mm
            
            if(self.t_as_i <= 0):
                print("ERROR: invalid t_as value")
                
            elif(self.t_as_i <= 0.1):
                print("Caution! Low value of t_as")
                
            #Vazão mássica ideal
            self.m_i = self.Mu_i*(np.pi/4*self.fr.d_n_2**2)*1e-6*np.sqrt(2*self.Delta_P*1e5*self.fr.Rho_2)*1e3 #vazão mássica ideal [g/s]
            
            #####   Cálculo do Estágio Real   #####
                        
            #Velocidade Tangencial
            self.W = self.m_i*1e-3/(self.fr.N_2*(np.pi/4*self.fr.d_in_2**2)*1e-6*self.fr.Rho_2) #Tangential Velocity

            
            
            #Número de Reynolds nos canais tangenciais
            self.Re_in = self.fr.Rho_2*self.W*self.fr.d_in_2*1e-3/self.fr.Din_visc_2
            
            

            if(self.Re_in < 2600): print("ERROR: Low Reynolds number (%.2E < 2600) [\u0394P = %.2f bar]"%(self.Re_in,self.Delta_P))
            elif(self.Re_in > 5e4): print("Caution: Reynolds number over 5e4")

            #Coeficiente de resistência de superfície
            self.Lambda = 0.3164*self.Re_in**(-0.25) #Bazarov
            #Lambda = 10**(25.8/(np.log10(Re_in))**2.58-2) #Bayvel
            
            
            
            #Parâmetro Geométrico equivalente
            self.A_w = self.calc_A (self.fr.R_0_2, (self.fr.d_n_2/2), self.fr.N_2, (self.fr.d_in_2/2), self.Lambda) 
            
            
            
            ##### Cálculo do Coeficiente de área livre real
            
        
            self.Phi_w = self.calc_Phi(self.A_w)
            
            #####
            
            #Coeficiente de Descarga
            
            self.Mu_w = self.Phi_w*np.sqrt(self.Phi_w)/(np.sqrt(2-self.Phi_w))
            
            
            #Espessura de camada líquida real
            self.t_as_w = self.fr.d_n_2/2*(1-np.sqrt(1-self.Phi_w)) #mm
            
            
            
            if(self.t_as_w <= 0):
                print("ERROR: invalid fluid layer thickness value")
                
            elif(self.t_as_w <= 0.1):
                print("Caution! Low fluid layer thickness")
            
            
            #Coeficiente de perda de momento angular
            self.K = self.A_w/self.A_i
            
            
            
            #Coeficiente de redução de área do estágio
            self.C = self.fr.R_0_2/(self.fr.d_n_2/2)
            
            
            ##### Cálculo da perda nos canais tangenciais #####
            
            
            self.Alpha_in = 90 - 180/np.pi*np.arctan(self.fr.d_s_2/2/self.fr.L_in_2)
            
            #Coeficiente de perda
            self.Ksi_in = -1*self.Alpha_in/150 + 1.1#-0.0204*(d_bx/d_in)**4-0.026*(d_bx/d_in)**3+0.3587*(d_bx/d_in)**2-0.6995*(d_bx/d_in)+0.5002
            
            
            
            #Fator de perda nos canais tangenciais
            self.Ksi = self.Ksi_in + self.Lambda*self.fr.L_in_2/(self.fr.d_in_2)
            
            
            
            #Vazão mássica real
            self.m_r = np.pi/4*(self.fr.d_n_2*1e-3)**2*np.sqrt(2*self.fr.Rho_2*self.Delta_P*1e5)/np.sqrt((2-self.Phi_w)/self.Phi_w**3+self.Ksi*self.A_i**2/self.fr.R_0_2**2)*1e3  #np.pi*(R_n*1e-3)**2*np.sqrt(2*Rho_p*Delta_P*1e5)/np.sqrt(1/Phi_w**2+A_i**2/(1-Phi_w)+(Ksi*N_in)*A_i**2/C**2)*1000
            
            #coeficiente de descarga total
            self.Mu_t = self.Mu_w/(np.sqrt(1+self.Mu_w**2*self.Ksi*self.A_i**2/self.C**2))
            
            
            self.Mass_ideal.append(self.m_i)
        
            
            
            self.Mass_real.append(self.m_r)
            
            self.Alpha2_i = 2*self.rad_to_deg * np.arcsin( 2*self.Mu_i*self.A_i / ((1+np.sqrt(1-self.Phi_i)) )
            
            self.Alpha2_w = 2*self.rad_to_deg * np.arcsin( 2*self.Mu_w*self.A_w / ((1+np.sqrt(1-self.Phi_i)) * np.sqrt(1-self.Ksi*self.Mu_w**2*self.A_i**2/self.C**2)) )
            
            
            
            #self.test.append(self.Alpha2_w)
            #self.test.append(2*180/np.pi*np.arcsin(2*self.Mu_i*self.A_w/(1+np.sqrt(1-self.Phi_i))*np.sqrt(1-self.Ksi*self.Mu_i**2*self.A_i**2/self.C**2)))
            #print(self.Alpha2_w - 2*180/np.pi*np.arcsin(2*self.Mu_i*self.A_w/(1+np.sqrt(1-self.Phi_i))*np.sqrt(1-self.Ksi*self.Mu_i**2*self.A_i**2/self.C**2)))
            self.test.append(2*self.Mu_w*self.A_w / ((1+np.sqrt(1-self.Phi_i)) * np.sqrt(1-self.Ksi*self.Mu_w**2*self.A_i**2/self.C**2)) )
            #self.test.append(2*self.Mu_i*self.A_w)
            #self.test.append((1+np.sqrt(1-self.Phi_i)))
            #self.test.append(np.sqrt(1-self.Ksi*self.Mu_i**2*self.A_i**2/self.C**2))
            
            #self.Alpha2_w = 2*180/np.pi*np.arcsin(2*self.Mu_t*self.A_i*self.K/(1+np.sqrt(1+self.Phi_w)*np.sqrt(1-self.Ksi*self.Mu_t**2*self.A_i**2)/self.C**2)) if -1 <= 2*self.Mu_t*self.A_i*self.K/(1+np.sqrt(1+self.Phi_w)*np.sqrt(1-self.Ksi*self.Mu_t**2*self.A_i**2)/self.C**2) <= 1 else 0
            
            self.Alpha2_ideal.append(self.Alpha2_i)
            
            self.Alpha2_real.append(self.Alpha2_w)
            
            self.Re_plot.append(self.Re_in)
            
            self.Ksi_plot.append(self.Ksi)

            
            self.M_target = []
            
            for self.i in self.dPressure:
                self.M_target.append(self.fr.m_2)
            
        for self.i in range(len(self.dPressure)):
            if(abs(self.dPressure[self.i] - self.fr.Delta_P0_2) < 0.005):

                print("\nm_bar_i = %.2f g/s"%(self.Mass_ideal[self.i]))
                print("\nm_bar = %.2f g/s"%(self.Mass_real[self.i]))
                print("%%target = %.2f %%"%(self.Mass_real[self.i]/self.fr.m_2*100))
            
            
        fig, ax21 = plt.subplots(figsize=(8, 6))
        ax21.grid()
        ax21.set_title("ST2 - ṁ x ΔP") 
        ax21.plot(self.dPressure, self.Mass_ideal, 'k--', markersize=1, label = "ṁ fluido ideal", color = '0.5')
        ax21.plot(self.dPressure, self.Mass_real, 'k:', markersize=3, label = "ṁ fluido viscoso")
        ax21.plot(self.dPressure, self.M_target, 'k-.', markersize=1, label = "ṁ alvo")
        ax21.legend()
        ax21.set_ylabel("ṁ [g/s]")
        ax21.set_xlabel("ΔP [bar]")
        ax21.set_xticks(np.arange(0,10.0+0.5,0.5))
        ax21.set_xlim(0.0,10.0)
        #ax21.grid()
        
        fig, ax22 = plt.subplots(figsize=(8, 6))
        ax22.grid()
        ax22.set_title("ST2 - Re x ΔP") 
        ax22.plot(self.dPressure, self.Re_plot, 'k-', markersize=1)
        ax22.set_ylabel("Re [ad.]")
        ax22.set_xlabel("ΔP [bar]")
        ax22.set_xticks(np.arange(0,10.0+0.5,0.5))
        ax22.set_xlim(0.0,10.0)
        
        fig, ax23 = plt.subplots(figsize=(8, 6))
        ax23.grid()
        ax23.set_title("ST2 - Ksi x ΔP") 
        ax23.plot(self.dPressure, self.Ksi_plot, 'k-', markersize=1)
        ax23.set_ylabel("Ksi [ad.]")
        ax23.set_xlabel("ΔP [bar]")
        ax23.set_xticks(np.arange(0,10.0+0.5,0.5))
        ax23.set_xlim(0.0,10.0)
            
        fig, ax24 = plt.subplots(figsize=(8, 6))
        ax24.grid()
        ax24.set_title("ST2 - 2α x ΔP") 
        ax24.plot(self.dPressure, self.Alpha2_ideal, 'k--', markersize=1, label = "2α fluido ideal", color = '0.5')
        ax24.plot(self.dPressure, self.Alpha2_real, 'k:', markersize=1, label = "2α fluido viscoso")
        ax24.legend()
        ax24.set_ylabel("2α [deg]")
        ax24.set_xlabel("ΔP [bar]")
        ax24.set_xticks(np.arange(0,10.0+0.5,0.5))
        ax24.set_xlim(0.0,10.0)
        
        fig, ax25 = plt.subplots(figsize=(8, 6))
        ax25.grid()
        ax25.set_title("ST2 - 2α x ΔP") 
        ax25.plot(self.dPressure, self.test, 'k--', markersize=1, label = "2α fluido ideal", color = '0.5')
        #ax25.plot(self.dPressure, self.Alpha2_real, 'k:', markersize=1, label = "2α fluido viscoso")
        ax25.legend()
        ax25.set_ylabel("2α [deg]")
        ax25.set_xlabel("ΔP [bar]")
        ax25.set_xticks(np.arange(0,10.0+0.5,0.5))
        ax25.set_xlim(0.0,10.0)
        