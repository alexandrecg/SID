# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 07:06:17 2020

@author: A. Goulart
"""

import numpy as np

import matplotlib.pyplot as plt

from M1_FileReader import FileReader_1





class Method_1:
    def __init__(self, foldername, filename):
        self.it_lim = 1000
        self.erro_max = 1e-3
        
        self.fr = FileReader_1()
        self.fr.setFolderName(foldername)
        self.fr.setFileName(filename)
        

        
    def set_it_lim(self,n):
        self.it_lim = n
        
    def get_it_lim(self):
        print("it_lim = ", self.it_lim)
        
    def set_erro_max(self,n):
        self.erro_max = n
        
    def get_erro_max(self):
        print("erro_max = ", self.erro_max)
    
    def run_M1(self):
        
        self.fr.read()
        
        print("Config: ", self.fr.Config)
        

        #if(self.fr.Config == "Mono\n"):
         #   recess = self.fr.recess
          #  t_w = self.fr.t_w
        
        ###### Method database generation ######
        

        self.Phi_min = 0.001
        
        self.Phi_max = 0.999
        
        self.Phi_step = (self.Phi_max-self.Phi_min)/self.it_lim
        
        self.Phi = np.arange(self.Phi_min,self.Phi_max+self.Phi_step,self.Phi_step)
        
        self.Mi = []
        
        self.A = []
        
        self.Alpha2 = []
        
        for phi in self.Phi:
            self.Mi_temp = phi*np.sqrt(phi/(2-phi))
            self.Alpha2_temp = 2*180/np.pi*np.arctan(np.sqrt(2*(1-phi)/phi))
            self.A_temp = (1-phi)/phi*np.sqrt(2/phi)
            
            self.Mi.append(self.Mi_temp)
            self.Alpha2.append(self.Alpha2_temp)
            self.A.append(self.A_temp)
            
        

        ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
        ##### ##### ##### Parte 1 - Cálculo do ST1 #### ##### ##### ##### ##### #####
        ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
        
        #self.Alpha2_1 = -1
        
        self.contador = 0
        
        self.Check = 0
        
        for i in range(self.it_lim):
            if(abs(self.Alpha2[i]-self.fr.Alpha2_1)<0.1):
                self.Mi_1 = self.Mi[i]
                self.A_1 = self.A[i]
                self.Phi_1 = self.Phi[i]
                self.Check = 1
        
        if(self.Check == 0): 
            print("ERROR: Mi_1")
            self.Check = 0
        
        self.R_n_1 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(self.fr.m_1*1e-3/(self.Mi_1*np.sqrt(2*self.fr.Rho_1*self.fr.Delta_P_1*1e5)))
        
        print("Initial Values: \n\tA_1 = \t%.1f; \n\tR_n_1 = \t%.3f \tmm; \n\t2Alpha = %.1f \tdeg"%(self.A_1,self.R_n_1,self.fr.Alpha2_1))

        
        self.entrance = []
        
        self.erro_R = 100
        
        while self.erro_R >= self.erro_max and self.contador < self.it_lim:
            

            self.R_in_1 = self.fr.Rs_Rn_1*self.R_n_1
            
            self.r_in_1 = np.sqrt(self.R_in_1*self.R_n_1/(self.fr.n_1*self.A_1))
            
            self.entrance.append(self.r_in_1)
            
            self.L_in_1 = self.fr.LD_in_1*self.r_in_1
            
            self.L_n_1 = self.fr.LD_n_1*self.R_n_1
            
            self.L_s_1 = self.fr.LD_s_1*self.R_in_1
            
            self.R_s_1 = self.R_in_1+self.r_in_1
        
            self.Re_in_1 = (2/np.pi)*self.fr.m_1*1e-3/(self.fr.n_1*self.r_in_1*1e-3*self.fr.Din_visc_1)
            
            self.Lambda_1 = 0.3164/(self.Re_in_1**0.25)
            
            self.Alpha_in_1 = 90 - 180/np.pi * np.arctan(self.R_s_1/self.L_in_1)
            
            
            if(self.fr.In_type_1 == "curved\n"):
                self.Ksi_in_1 = -1*self.Alpha_in_1/150 + 1.1
                
            elif(self.fr.In_type_1 == "straigth\n"):
                self.Ksi_in_1 = 0.5*np.exp(-1.4*(self.fr.entrance_radius_1/(2*r_in_1)))
                
            else:
                print("ERROR: Injector inlet type not recognized on Stage 1")
            
            self.Ksi_1 = self.Ksi_in_1 + self.Lambda_1*self.L_in_1/(2*self.r_in_1)
            
            
            self.A_eq_1 = self.R_in_1*self.R_n_1/(self.fr.n_1*self.r_in_1**2+self.Lambda_1/2*self.R_in_1*(self.R_in_1-self.R_n_1))
            
            for i in range(self.it_lim):
                if(abs(self.A[i]-self.A_eq_1)<0.05):
                    self.Mi_eq_1 = self.Mi[i]
                    self.Alpha2_eq_1 = self.Alpha2[i]
                    self.Phi_eq_1 = self.Phi[i]
                    self.Check = 1
            
            if(self.Check == 0): 
                print("ERROR: Mi_eq_1")
                self.Check = 0
            
            
            self.Mi_i_1 = self.Mi_eq_1/np.sqrt(1+self.Ksi_1*self.Mi_eq_1**2*self.A_1**2/self.fr.Rs_Rn_1**2)
            
            self.Alpha2_eq_1_calc = 2*180/np.pi*np.arcsin(2*self.Mi_i_1*self.A_eq_1/((1+np.sqrt(1-self.Phi_eq_1))*np.sqrt(1-self.Ksi_1*self.Mi_i_1**2*self.A_1**2/self.fr.Rs_Rn_1**2)))
            
            self.erro_R = self.R_n_1
            
            self.R_n_1 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(self.fr.m_1*1e-3/(self.Mi_i_1*np.sqrt(2*self.fr.Rho_1*self.fr.Delta_P_1*1e5)))
            
            self.erro_R = abs(self.erro_R-self.R_n_1)
            
            self.A_1 = self.R_in_1*self.R_n_1/(self.fr.n_1*self.r_in_1**2)
            
            self.r_mn_1 = self.R_n_1*np.sqrt(1-self.Phi_eq_1)
            
            self.t_fluid_1 = self.R_n_1-self.r_mn_1
            
            self.r_mk_1 = self.r_mn_1 * np.sqrt(2*(1-self.Phi_eq_1)/(2-self.Phi_eq_1))
            
            
            print("\n\nIt. %d) \n\tA_1 = \t%.2f; \n\tR_n_1 = \t%.3f \tmm; \n\t2Alpha = %.1f (%.1f) \tdeg"%(self.contador+1,self.A_1,self.R_n_1,self.Alpha2_eq_1,self.Alpha2_eq_1_calc))
        
            self.contador += 1
        
        fig, ax1 = plt.subplots()
        ax1.plot(range(self.contador),self.entrance)
        ax1.set_title("Stage 1 - Solution convergence (r_in)")
        ax1.set_xlabel("Iteration")
        ax1.set_ylabel("r_in [mm]")
        ax1.grid()
        
        print("\nPhi = %.2f;Phi_eq = %.2f \n"%(self.Phi_1,self.Phi_eq_1))
        
        print("r_mn_1 = %.2f mm"%(self.r_mn_1))
        print("t_fluid_1 = %.2f mm"%(self.t_fluid_1))
        
        print("\nr_mk_1 = %.2f mm"%(self.r_mk_1))
        
        
        print("\n\n\t ST1 Injector Geometry: \t")
        print("Number of Inlet Channels = \t\t%.1f"%(self.fr.n_1))
        
        print("\nRadial Dimensions:")
        print("Nozzle Radius = \t\t\t%.3f \tmm"%(self.R_n_1))
        print("Inlet Radius =  \t\t\t%.3f \tmm"%(self.R_in_1))
        print("Swirl chamber Radius = \t\t%.3f \tmm"%(self.R_s_1))
        print("Inlet Channel Radius = \t\t%.3f \tmm"%(self.r_in_1))
        
        print("\nLinear dimensions:")
        print("Nozzle Length = \t\t\t%.3f \tmm"%(self.L_n_1))
        print("Swirl Chamber Length = \t\t%.3f \tmm"%(self.L_s_1))
        print("Inlet Channel Length = \t\t%.3f \tmm"%(self.L_in_1))
        
        print("\n\nMass flow target = \t\t%.2f g/s"%(self.fr.m_1))
        print("Delta P ST1 = \t\t%.2f bar\n"%(self.fr.Delta_P_1))
        
        if(np.pi*self.R_n_1**2 < self.fr.n_1*np.pi*self.r_in_1**2):
            print("Warning: Exit area smaller than Entrance area.")
        
        #Cálculo da Não-uniformidade
        
        I_1 = 23.7/((self.R_in_1/self.R_n_1)**2.7*self.fr.n_1**1.34*self.Phi_eq_1**1.1*(self.L_s_1/(2*self.R_s_1))**0.15)
        
        print("\n\nNão Uniformidade Esperada = \t%.2f %%"%(I_1))
        
        
        ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
        ##### ##### ##### Parte 2 - Cálculo do ST2 #### ##### ##### ##### ##### #####
        ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
        
        if(self.fr.Config == "Bi\n"):
            
            self.contador_max = 10
            
            self.contador = 0
            
            self.Check = 0
            
            for i in range(self.it_lim):
                if(abs(self.Alpha2[i]-self.fr.Alpha2_2)<0.1):
                    self.Mi_2 = self.Mi[i]
                    self.A_2 = self.A[i]
                    self.Phi_2 = self.Phi[i]
                    self.Check = 1
            
            if(self.Check == 0): 
                print("ERROR: Mi_2")
                self.Check = 0
            
            self.R_n_2 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(self.fr.m_2*1e-3/(self.Mi_2*np.sqrt(2*self.fr.Rho_2*self.fr.Delta_P_2*1e5)))
            
            print("\nA_2 = \t%.1f; R_n_2 = \t%.3f \tmm; 2Alpha = %.1f \tdeg"%(self.A_2,self.R_n_2,self.fr.Alpha2_2))
            
            self.entrance = []
            
            self.R_in_2 = self.fr.Rs_Rn_2*self.R_n_2
            
            self.r_in_2_0 = np.sqrt(self.R_in_2*self.R_n_2/(self.fr.n_2*self.A_2))
            
            self.r_in_2 = self.r_in_2_0
            
            
            while self.contador < self.contador_max:
                    
                self.fr.Rs_Rn_2
                
                #self.R_in_2 = self.R_n_2-self.r_in_2
                self.R_in_2 = self.fr.Rs_Rn_2*self.R_n_2
                
                self.r_in_2 = np.sqrt(self.R_in_2*self.R_n_2/(self.fr.n_2*self.A_2))
                
                self.entrance.append(self.r_in_2)
                
                self.L_in_2 = self.fr.LD_in_2*self.r_in_2
                
                self.L_n_2 = self.fr.LD_n_2*self.R_n_2
                
                self.L_s_2 = self.fr.LD_s_2*self.R_in_2
                
                #self.R_s_2 = self.R_n_2 
                self.R_s_2 = self.R_in_2+self.r_in_2
            
                self.Re_in_2 = (2/np.pi)*self.fr.m_2*1e-3/(self.fr.n_2*self.r_in_2*1e-3*self.fr.Din_visc_2)
                
                self.Lambda_2 = 0.3164/(self.Re_in_2**0.25)
                
                self.Alpha_in_2 = 90 - 180/np.pi * np.arctan(self.R_s_2/self.L_in_2)
                
                
                if(self.fr.In_type_2 == "curved\n"):
                    self.Ksi_in_2 = -1*self.Alpha_in_2/150 + 1.1
                    
                elif(self.fr.In_type_2 == "straigth\n"):
                    self.Ksi_in_2 = 0.5*np.exp(-1.4*(self.fr.entrance_radius_2/(2*r_in_2)))
                    
                else:
                    print("ERROR: Injector inlet type not recognized on Stage 2")
                
                self.Ksi_2 = self.Ksi_in_2 + self.Lambda_2*self.L_in_2/(2*self.r_in_2)
                
                
                self.A_eq_2 = self.R_in_2*self.R_n_2/(self.fr.n_2*self.r_in_2**2+self.Lambda_2/2*self.R_in_2*(self.R_in_2-self.R_n_2))
                
                for i in range(self.it_lim):
                    if(abs(self.A[i]-self.A_eq_2)<0.05):
                        self.Mi_eq_2 = self.Mi[i]
                        self.Alpha2_eq_2 = self.Alpha2[i]
                        self.Phi_eq_2 = self.Phi[i]
                        self.Check = 1
                
                if(self.Check == 0): 
                    print("ERROR: Mi_eq_2")
                    self.Check = 0
                
                #R_fluid = 
                
                self.Mi_i_2 = self.Mi_eq_2/np.sqrt(1+self.Ksi_2*self.Mi_eq_2**2*self.A_2**2/self.fr.Rs_Rn_2**2)
                
                self.Alpha2_eq_2_calc = 2*180/np.pi*np.arcsin(2*self.Mi_i_2*self.A_eq_2/((1+np.sqrt(1-self.Phi_eq_2))*np.sqrt(1-self.Ksi_2*self.Mi_i_2**2*self.A_2**2/self.fr.Rs_Rn_2**2)))
                
                self.R_n_2 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(self.fr.m_2*1e-3/(self.Mi_i_2*np.sqrt(2*self.fr.Rho_2*self.fr.Delta_P_2*1e5)))
                
                self.A_2 = self.R_in_2*self.R_n_2/(self.fr.n_2*self.r_in_2**2)
                
                self.r_mn_2 = self.R_n_2*np.sqrt(1-self.Phi_eq_2)
                
                self.t_fluid_2 = self.R_n_2-self.r_mn_2
                
                self.r_mk_2 = self.r_mn_2 * np.sqrt(2*(1-self.Phi_eq_2)/(2-self.Phi_eq_2)) #R_n_2*np.sqrt(2*(1-Phi_eq_2)**2/(2-Phi_eq_2))
                
                
                print("\n\nIt. %d) A_2 = \t%.2f; R_n_2 = \t%.3f \tmm; 2Alpha = %.1f (%.1f) \tdeg"%(self.contador+1,self.A_2,self.R_n_2,self.Alpha2_eq_2,self.Alpha2_eq_2_calc))
            
                self.contador += 1
            
            fig, ax2 = plt.subplots()
            ax2.plot(range(self.contador_max),self.entrance)
            ax2.set_title("Stage 2 - Solution convergence (r_in)")
            ax2.set_xlabel("Iteration")
            ax2.set_ylabel("r_in [mm]")
            ax2.grid()
            
            print("\nPhi = %.2f;Phi_eq = %.2f \n"%(self.Phi_2,self.Phi_eq_2))
            
            print("r_mn_2 = %.2f mm"%(self.r_mn_2))
            print("t_fluid_2 = %.2f mm"%(self.t_fluid_2))
            
            print("\nr_mk_2 = %.2f mm"%(self.r_mk_2))
            
            print("\n\n\t ST2 Injector Geometry: \t")
            print("Number of Inlet Channels = \t\t%.1f"%(self.fr.n_2))
            
            print("\nRadial Dimensions:")
            print("Nozzle Radius = \t\t\t%.3f \tmm"%(self.R_n_2))
            print("Inlet Radius =  \t\t\t%.3f \tmm"%(self.R_in_2))
            print("Swirl chamber Radius = \t\t%.3f \tmm"%(self.R_s_2))
            print("Inlet Channel Radius = \t\t%.3f \tmm"%(self.r_in_2))
            
            if(self.R_n_2 > self.R_s_2):
                print("\n Erro de dimensão: R_n_2 > R_s_2 não é permitido. Favor aumentar o valor de Rs_Rn_2.\n")
            
            print("\nLinear dimensions:")
            print("Nozzle Length = \t\t\t%.3f \tmm"%(self.L_n_2))
            print("Swirl Chamber Length = \t\t%.3f \tmm"%(self.L_s_2))
            print("Inlet Channel Length = \t\t%.3f \tmm"%(self.L_in_2))
            
            print("\n\nMass flow target = \t\t%.2f g/s"%(self.fr.m_2))
            print("Delta P ST2 = \t\t%.2f bar\n"%(self.fr.Delta_P_2))
            
            
            #Cálculo da Não-uniformidade
        
            I_2 = 23.7/((self.R_in_2/self.R_n_2)**2.7*self.fr.n_2**1.34*self.Phi_eq_2**1.1*(self.L_s_2/(2*self.R_s_2))**0.15)
        
            print("\n\nNão Uniformidade Esperada = \t%.2f %%"%(I_2))
        
        
            ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
            ##### ##### ##### Parte 3 - Integração dos Estágios## ##### ##### ##### #####
            ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
            
            print("\n\n\t Geometry Check: \t")
            
            print("D_n_1 = %.3f mm"%(2*self.R_n_1))
            
            print("D_w_1 = %.3f mm"%(2*self.R_n_1+2*self.fr.t_w))
            
            print("D_n_2 = %.3f mm"%(2*self.R_n_2))
            
            self.Separation = self.r_mk_2 - (self.R_n_1 + self.fr.t_w)
            
            self.fluid_t_2_s = self.R_s_2 - self.r_mk_2
            
            print("fluid_t_2_s = %.3f mm"%(self.fluid_t_2_s))
            
            #if(self.R_n_2 - self.R_n_1 - self.fr.t_w >= (1+self.fr.delta/100)*self.t_fluid_2):
            if(self.Separation >= self.fr.delta/100*self.fluid_t_2_s):
                print("GO! DeltaR = %.3f mm (mín: %.2f mm)"%(self.Separation, self.fr.delta/100*self.fluid_t_2_s))
                
            else:
                print("NOGO! DeltaR = %.3f mm (mín: %.2f mm)"%(self.Separation, self.fr.delta/100*self.fluid_t_2_s))
            
            ##### #####
            # Alteração do ângulo do ST2 #
            ##### #####
            
            print("\n\n\t ST2 Modified Injector Geometry: \t")
            
            self.Alpha2_2_goal = self.Alpha2_eq_1_calc - self.fr.angle_dif
            
            print("2Alpha_2_goal = \t%.1f deg"%(self.Alpha2_2_goal))
            
            self.a_2 = 2*(1-self.Phi_eq_2)**2/(2-self.Phi_eq_2)
            
            self.R_out_2 = self.R_n_2*np.sqrt(self.a_2*(1+(np.tan(self.Alpha2_2_goal/2*np.pi/180))**-2))
            
            print("R_out_2 = \t\t\t%.3f mm"%(self.R_out_2))
            
            self.L_out_2 = self.fr.LD_n_3*self.R_n_2#(self.R_out_2-self.R_n_1-self.fr.t_w)
            
            print("L_out_2 = \t\t\t%.3f"%(self.L_out_2))
            
            self.Phi_out = self.Phi_eq_2/(3-2*self.Phi_eq_2)
            
            self.R_out_mn = self.R_out_2*np.sqrt(1-self.Phi_out)
            
            print("R_out_mn = \t\t\t%.3f mm"%(self.R_out_mn))
            
            self.t_fluid_out = self.R_out_2-self.R_out_mn
            
            print("t_fluid_out = \t\t%.3f mm"%(self.t_fluid_out))
            
            self.inj_3_H = [0,self.L_out_2]
            self.inj_3_D_in = [self.R_out_2,self.R_out_2]
            self.inj_3_D_in_mirror = [-1*self.R_out_2,-1*self.R_out_2]
            
            self.y_i_3 = 0
            self.y_f_3 = -5
            self.y_step_3 = (self.y_f_3-self.y_i_3)/100
            self.y_3 = np.arange(self.y_i_3,self.y_f_3+self.y_step_3,self.y_step_3)
            
            self.x_3 = []
            self.x_3_mirror = []
            
            self.R_av_3 = (self.R_out_2 + self. R_out_mn)/2
            
            for y in self.y_3:
                self.x_temp = self.R_av_3*np.sqrt(1+(y)**2/(self.R_av_3*np.tan((90-0.5*self.Alpha2_2_goal)*np.pi/180))**2)
                self.x_3.append(self.x_temp)
                self.x_3_mirror.append(-1*self.x_temp)
            
            self.H_1_n_top = self.fr.recess+self.L_n_1
                
            self.H_1_trans = (self.R_s_1-self.R_n_1)/np.tan(self.fr.Trans_angle_1*np.pi/180)
            
            self.H_1_s_top = self.H_1_n_top + self.L_s_1 - self.H_1_trans
            
            self.H_2_n_top = self.L_n_2
    
            self.H_2_trans = (self.R_s_2-self.R_n_2)/np.tan(self.fr.Trans_angle_2*np.pi/180)
        
            self.H_2_s_top = self.H_2_n_top + self.L_s_2 - self.H_2_trans
            
            self.inj_D_top = [-1*self.R_s_1,self.R_s_1]
            self.inj_H_top = [self.H_1_s_top,self.H_1_s_top]
            
            
            
            ##### #####
            
            self.inj_H = [self.fr.recess,self.L_n_1+self.fr.recess]
            self.inj_D_in = [self.R_n_1,self.R_n_1]
            self.inj_D_in_mirror = [-1*self.R_n_1,-1*self.R_n_1]
            
            
            self.inj_H_out = [self.fr.recess,self.H_2_s_top]
            self.inj_D_out = [self.R_n_1+self.fr.t_w,self.R_n_1+self.fr.t_w]
            self.inj_D_out_mirror = [-1*self.R_n_1-self.fr.t_w,-1*self.R_n_1-self.fr.t_w]
            
            #### #### ####
    
            
            self.inj_1_Strans_D = [self.R_n_1,self.R_s_1]
            self.inj_1_Strans_D_mirror = [-1*self.R_n_1,-1*self.R_s_1]
            self.inj_1_Strans_H = [self.H_1_n_top,self.H_1_n_top + self.H_1_trans]
            
            self.inj_1_Sfluid_D = [self.r_mk_1,self.r_mk_1]
            self.inj_1_Sfluid_D_mirror = [-1*self.r_mk_1,-1*self.r_mk_1]
            
            self.inj_1_S_D = [self.R_s_1,self.R_s_1]
            self.inj_1_S_D_mirror = [-1*self.R_s_1,-1*self.R_s_1]
            self.inj_1_S_H = [self.H_1_n_top + self.H_1_trans,self.H_1_s_top]
            
            self.inj_1_Sfluid_Dtrans = [self.r_mn_1,self.r_mk_1]
            self.inj_1_Sfluid_Dtrans_mirror = [-1*self.r_mn_1,-1*self.r_mk_1]
            self.inj_1_Sfluid_Htrans = [self.H_1_n_top,self.H_1_n_top + self.H_1_trans]
            
            self.inj_1_in_top = []
            self.inj_1_in_bottom = []
            
            self.circ_x_1 = np.arange(-1*self.R_s_1,-1*self.R_s_1+2*self.r_in_1+0.005,0.005)
            
            for x in self.circ_x_1:
                self.inj_1_in_top.append((self.H_1_s_top-self.r_in_1) + np.sqrt(self.r_in_1**2-(x+self.R_s_1-self.r_in_1)**2))
                self.inj_1_in_bottom.append((self.H_1_s_top-self.r_in_1) - np.sqrt(self.r_in_1**2-(x+self.R_s_1-self.r_in_1)**2))
                    
            #### #### ####
            
            self.inj_F_H = [0+self.fr.recess,0+self.fr.recess]
            self.inj_F = [self.R_n_1,self.R_n_1+self.fr.t_w]
            self.inj_F_mirror = [-1*self.R_n_1,-1*self.R_n_1-self.fr.t_w]
            
            if(self.fr.Config == "Bi\n"):
                
                
                
                self.inj_D_2_top = [self.R_s_2,(self.R_n_1+self.fr.t_w)]
                self.inj_D_2_top_mirror = [-1*self.R_s_2,-1*(self.R_n_1+self.fr.t_w)]
                self.inj_H_2_top = [self.H_2_s_top,self.H_2_s_top]
                
                
                self.inj_2_H = [0,self.L_n_2]
                self.inj_2_D_in = [self.R_n_2,self.R_n_2]
                self.inj_2_D_in_mirror = [-1*self.R_n_2,-1*self.R_n_2]
                self.inj_2_D_out = [self.R_n_2+2*self.fr.t_w,self.R_n_2+2*self.fr.t_w]
                self.inj_2_D_out_mirror = [-1*self.R_n_2-2*self.fr.t_w,-1*self.R_n_2-2*self.fr.t_w]
                
                self.inj_2_F_H = [0,0]
                self.inj_2_F = [self.R_n_2,self.R_n_2+2*self.fr.t_w]
                self.inj_2_F_mirror = [-1*self.R_n_2,-1*self.R_n_2-2*self.fr.t_w]
            
            
            
            
            self.fluid_D_1 = [self.r_mn_1,self.r_mn_1]
            self.fluid_D_1_mirror = [-1*self.r_mn_1,-1*self.r_mn_1]
            
            self.y_i_1 = 0+self.fr.recess
            self.y_f_1 = -5
            self.y_step_1 = (self.y_f_1-self.y_i_1)/100
            self.y_1 = np.arange(self.y_i_1,self.y_f_1+self.y_step_1,self.y_step_1)
            
            self.x_1 = []
            self.x_1_mirror = []
            
            self.R_av_1 = (self.R_n_1 + self.r_mn_1)/2
            
            for y in self.y_1:
                self.x_temp = self.R_av_1*np.sqrt(1+(y-self.fr.recess)**2/(self.R_av_1*np.tan((90-0.5*self.Alpha2_eq_1_calc)*np.pi/180))**2)
                self.x_1.append(self.x_temp)
                self.x_1_mirror.append(-1*self.x_temp)
            
            if(self.fr.Config == "Bi\n"):
                self.fluid_D_2 = [self.r_mn_2,self.r_mn_2]
                self.fluid_D_2_mirror = [-1*self.r_mn_2,-1*self.r_mn_2]
                
                self.y_i_2 = 0
                self.y_f_2 = -5
                self.y_step_2 = (self.y_f_2-self.y_i_2)/100
                self.y_2 = np.arange(self.y_i_2,self.y_f_2+self.y_step_2,self.y_step_2)
                
                self.x_2 = []
                self.x_2_mirror = []
                
                self.R_av_2 = (self.R_n_2 + self.r_mn_2)/2
    
                
                print("\n")
                print("1: ",self.Alpha2_eq_1_calc)
                print("2: ",self.Alpha2_eq_2_calc)
                print("3: ", self.Alpha2_2_goal)
            
            
                for y in self.y_2:
                    self.x_temp = self.R_av_2*np.sqrt(1+y**2/(self.R_av_2*np.tan((90-0.5*self.Alpha2_eq_2_calc)*np.pi/180))**2)
                    self.x_2.append(self.x_temp)
                    self.x_2_mirror.append(-1*self.x_temp)
            
            
            
                self.inj_2_Strans_D = [self.R_n_2,self.R_s_2]
                self.inj_2_Strans_D_mirror = [-1*self.R_n_2,-1*self.R_s_2]
                self.inj_2_Strans_H = [self.H_2_n_top,self.H_2_n_top + self.H_2_trans]
                
                self.inj_2_Sfluid_D = [self.r_mk_2,self.r_mk_2]
                self.inj_2_Sfluid_D_mirror = [-1*self.r_mk_2,-1*self.r_mk_2]
                
                self.inj_2_S_D = [self.R_s_2,self.R_s_2]
                self.inj_2_S_D_mirror = [-1*self.R_s_2,-1*self.R_s_2]
                self.inj_2_S_H = [self.H_2_n_top + self.H_2_trans,self.H_2_s_top]
                
                self.inj_2_Sfluid_Dtrans = [self.r_mn_2,self.r_mk_2]
                self.inj_2_Sfluid_Dtrans_mirror = [-1*self.r_mn_2,-1*self.r_mk_2]
                self.inj_2_Sfluid_Htrans = [self.H_2_n_top,self.H_2_n_top + self.H_2_trans]
                
                self.inj_2_in_top = []
                self.inj_2_in_bottom = []
                
                self.circ_x_2 = np.arange(-1*self.R_s_2,-1*self.R_s_2+2*self.r_in_2+0.005,0.005)
                
                for x in self.circ_x_2:
                    self.inj_2_in_top.append((self.H_2_s_top-self.r_in_2) + np.sqrt(self.r_in_2**2-(x+self.R_s_2-self.r_in_2)**2))
                    self.inj_2_in_bottom.append((self.H_2_s_top-self.r_in_2) - np.sqrt(self.r_in_2**2-(x+self.R_s_2-self.r_in_2)**2))
            
            
            
            fig, ax3 = plt.subplots()
            ax3.grid()
    
            
            ax3.plot(self.inj_D_in,self.inj_H, linewidth = 2, color = 'red')
            ax3.plot(self.inj_D_in_mirror,self.inj_H, linewidth = 2, color = 'red')
            ax3.plot(self.inj_D_out,self.inj_H_out, linewidth = 2, color = 'red')
            ax3.plot(self.inj_D_out_mirror,self.inj_H_out, linewidth = 2, color = 'red')
            ax3.plot(self.inj_F,self.inj_F_H, linewidth = 2, color = 'red')
            ax3.plot(self.inj_F_mirror,self.inj_F_H, linewidth = 2, color = 'red')
            
            ax3.plot(self.inj_D_top,self.inj_H_top, linewidth = 2, color = 'red')
            
            ax3.plot(self.inj_1_Sfluid_D,self.inj_1_S_H, linewidth = 1, color = 'red', linestyle = '--')
            ax3.plot(self.inj_1_Sfluid_D_mirror,self.inj_1_S_H, linewidth = 1, color = 'red', linestyle = '--')
            ax3.plot(self.inj_1_Sfluid_Dtrans,self.inj_1_Sfluid_Htrans, linewidth = 1, color = 'red', linestyle = '--')
            ax3.plot(self.inj_1_Sfluid_Dtrans_mirror,self.inj_1_Sfluid_Htrans, linewidth = 1, color = 'red', linestyle = '--')
            
            
            ax3.plot(self.inj_1_S_D,self.inj_1_S_H, linewidth = 2, color = 'red')
            ax3.plot(self.inj_1_S_D_mirror,self.inj_1_S_H, linewidth = 2, color = 'red')
            ax3.plot(self.inj_1_Strans_D,self.inj_1_Strans_H, linewidth = 2, color = 'red')
            ax3.plot(self.inj_1_Strans_D_mirror,self.inj_1_Strans_H, linewidth = 2, color = 'red')
            
            ax3.plot(self.circ_x_1,self.inj_1_in_top, linewidth = 1, color = 'red')
            ax3.plot(self.circ_x_1,self.inj_1_in_bottom, linewidth = 1, color = 'red')
            
    
            ax3.plot(self.fluid_D_1,self.inj_H, linewidth = 1, color = 'red', linestyle = '--')
            ax3.plot(self.fluid_D_1_mirror,self.inj_H, linewidth = 1, color = 'red', linestyle = '--')
            ax3.plot(self.x_1,self.y_1, linewidth = 1, color = 'red', linestyle = '--')
            ax3.plot(self.x_1_mirror,self.y_1, linewidth = 1, color = 'red', linestyle = '--')
            
            ax3.plot(self.inj_D_2_top,self.inj_H_2_top, linewidth = 2, color = 'blue')
            ax3.plot(self.inj_D_2_top_mirror,self.inj_H_2_top, linewidth = 2, color = 'blue')
            
            
            ax3.plot(self.inj_2_Sfluid_D,self.inj_2_S_H, linewidth = 1, color = 'blue', linestyle = '--')
            ax3.plot(self.inj_2_Sfluid_D_mirror,self.inj_2_S_H, linewidth = 1, color = 'blue', linestyle = '--')
            ax3.plot(self.inj_2_Sfluid_Dtrans,self.inj_2_Sfluid_Htrans, linewidth = 1, color = 'blue', linestyle = '--')
            ax3.plot(self.inj_2_Sfluid_Dtrans_mirror,self.inj_2_Sfluid_Htrans, linewidth = 1, color = 'blue', linestyle = '--')
            
            
            ax3.plot(self.inj_2_S_D,self.inj_2_S_H, linewidth = 2, color = 'blue')
            ax3.plot(self.inj_2_S_D_mirror,self.inj_2_S_H, linewidth = 2, color = 'blue')
            ax3.plot(self.inj_2_Strans_D,self.inj_2_Strans_H, linewidth = 2, color = 'blue')
            ax3.plot(self.inj_2_Strans_D_mirror,self.inj_2_Strans_H, linewidth = 2, color = 'blue')
            
            ax3.plot(self.circ_x_2,self.inj_2_in_top, linewidth = 1, color = 'blue')
            ax3.plot(self.circ_x_2,self.inj_2_in_bottom, linewidth = 1, color = 'blue')
            
            ax3.plot(self.inj_2_D_in,self.inj_2_H, linewidth = 2, color = 'blue')
            ax3.plot(self.inj_2_D_in_mirror,self.inj_2_H, linewidth = 2, color = 'blue')
            ax3.plot(self.inj_2_F,self.inj_2_F_H, linewidth = 2, color = 'blue')
            ax3.plot(self.inj_2_F_mirror,self.inj_2_F_H, linewidth = 2, color = 'blue')
            
            ax3.plot(self.fluid_D_2,self.inj_2_H, linewidth = 1, color = 'blue', linestyle = '--')
            ax3.plot(self.fluid_D_2_mirror,self.inj_2_H, linewidth = 1, color = 'blue', linestyle = '--')
    
            ax3.plot(self.fluid_D_2,self.inj_2_H, linewidth = 1, color = 'blue', linestyle = '--')
            ax3.plot(self.fluid_D_2_mirror,self.inj_2_H, linewidth = 1, color = 'blue', linestyle = '--')
            ax3.plot(self.x_2,self.y_2, linewidth = 1, color = 'blue', linestyle = '--')
            ax3.plot(self.x_2_mirror,self.y_2, linewidth = 1, color = 'blue', linestyle = '--')
            
            #ax3.plot(inj_3_D_in,inj_3_H, linewidth = 2, color = '0.3')
            #ax3.plot(inj_3_D_in_mirror,inj_3_H, linewidth = 2, color = '0.3')
            
            #ax3.plot(x_3,y_3, linewidth = 1, color = '0.3', linestyle = '--')
            #ax3.plot(x_3_mirror,y_3, linewidth = 1, color = '0.3', linestyle = '--')
            
            ax3.set_xlabel("[mm]")
            ax3.set_ylabel("[mm]")
            
            self.min_x = min([-1*(self.R_n_2+self.fr.t_w+2.0),-1*(self.R_s_1+2.0)])
            self.max_x = max([(self.R_n_2+self.fr.t_w+2.0),(self.R_s_1+2.0)])
            self.min_y = -5.0
            self.max_y = self.H_1_s_top+5.0
            
            #ax3.axis(xlim=(self.min_x, self.max_x), ylim=(self.min_y, self.max_y))
            ax3.axis('equal')
            
            
            ##### ##### ##### ##### #####
            ##### ##### ##### ##### #####
            
            #### mod ####
            if(self.fr.Config == "Bi\n"):
                self.trans = 0.5
                
                self.inj_2_H = [self.L_out_2,self.L_out_2+self.L_n_2]
                self.inj_2_D_in = [self.R_n_2,self.R_n_2]
                self.inj_2_D_in_mirror = [-1*self.R_n_2,-1*self.R_n_2]
                #self.inj_2_D_out = [self.R_n_2+2*self.fr.t_w,self.R_n_2+2*self.fr.t_w]
                #self.inj_2_D_out_mirror = [-1*self.R_n_2-2*self.fr.t_w,-1*self.R_n_2-2*self.fr.t_w]
                
                self.inj_2_F_H = [0,0]
                self.inj_2_F = [self.R_out_2,self.R_out_2+2*self.fr.t_w]
                self.inj_2_F_mirror = [-1*self.R_out_2,-1*self.R_out_2-2*self.fr.t_w]
                
                ###
    
                self.H_2_n_top = self.L_out_2+self.L_n_2
            
                self.H_2_trans = (self.R_s_2-self.R_n_2)/np.tan(self.fr.Trans_angle_2*np.pi/180)
            
                self.H_2_s_top = self.H_2_n_top + self.L_s_2 - self.H_2_trans
    
        
                self.inj_H_out = [self.fr.recess,self.H_2_s_top]
                self.inj_D_out = [self.R_n_1+self.fr.t_w,self.R_n_1+self.fr.t_w]
                self.inj_D_out_mirror = [-1*self.R_n_1-self.fr.t_w,-1*self.R_n_1-self.fr.t_w]
    
    
                self.inj_D_2_top = [self.R_s_2,(self.R_n_1+self.fr.t_w)]
                self.inj_D_2_top_mirror = [-1*self.R_s_2,-1*(self.R_n_1+self.fr.t_w)]
                self.inj_H_2_top = [self.H_2_s_top,self.H_2_s_top]
    
                self.inj_2_Strans_D = [self.R_n_2,self.R_s_2]
                self.inj_2_Strans_D_mirror = [-1*self.R_n_2,-1*self.R_s_2]
                self.inj_2_Strans_H = [self.H_2_n_top,self.H_2_n_top + self.H_2_trans]
                
                self.inj_2_Sfluid_D = [self.r_mk_2,self.r_mk_2]
                self.inj_2_Sfluid_D_mirror = [-1*self.r_mk_2,-1*self.r_mk_2]
                
                self.inj_2_S_D = [self.R_s_2,self.R_s_2]
                self.inj_2_S_D_mirror = [-1*self.R_s_2,-1*self.R_s_2]
                self.inj_2_S_H = [self.H_2_n_top + self.H_2_trans,self.H_2_s_top]
                
                self.inj_2_Sfluid_Dtrans = [self.r_mn_2,self.r_mk_2]
                self.inj_2_Sfluid_Dtrans_mirror = [-1*self.r_mn_2,-1*self.r_mk_2]
                self.inj_2_Sfluid_Htrans = [self.H_2_n_top,self.H_2_n_top + self.H_2_trans]
                
                self.inj_2_in_top = []
                self.inj_2_in_bottom = []
                
                self.circ_x_2 = np.arange(-1*self.R_s_2,-1*self.R_s_2+2*self.r_in_2+0.005,0.005)
                
                for x in self.circ_x_2:
                    self.inj_2_in_top.append((self.H_2_s_top-self.r_in_2) + np.sqrt(self.r_in_2**2-(x+self.R_s_2-self.r_in_2)**2))
                    self.inj_2_in_bottom.append((self.H_2_s_top-self.r_in_2) - np.sqrt(self.r_in_2**2-(x+self.R_s_2-self.r_in_2)**2))
            
            
                
                self.H_3_trans = (self.R_out_2-self.R_n_2)/np.tan(self.fr.Trans_angle_3*np.pi/180)
                
                self.inj_3_Dtrans_in = [self.R_out_2,self.R_n_2]
                self.inj_3_Dtrans_in_mirror = [-1*self.R_out_2,-1*self.R_n_2]
                self.inj_3_Htrans = [self.L_out_2-self.H_3_trans,self.L_out_2]
                
                ###
                
                self.inj_3_H = [0,self.L_out_2-self.H_3_trans]
                self.inj_3_D_in = [self.R_out_2,self.R_out_2]
                self.inj_3_D_in_mirror = [-1*self.R_out_2,-1*self.R_out_2]
                
                
                #self.inj_23_H = [self.L_out_2-self.trans,self.L_out_2+self.trans]
                #self.inj_23_D_in = [self.R_out_2,self.R_n_2]
                #self.inj_23_D_in_mirror = [-1*self.R_out_2,-1*self.R_n_2]
                
                self.fluid_D_3 = [self.R_out_mn,self.R_out_mn]
                self.fluid_D_3_mirror = [-1*self.R_out_mn,-1*self.R_out_mn]
                
                self.fluid_trans = [self.r_mn_2,self.R_out_mn]
                self.fluid_trans_mirror = [-1*self.r_mn_2,-1*self.R_out_mn]
                self.fluid_trans_H = [self.L_out_2,self.L_out_2-self.H_3_trans]
                
                #for self.r_trans in np.arange(self.r_mn_2,self.R_out_mn,abs(self.r_mn_2-self.R_out_mn)/100):
                    #self.Phi_trans = self.Phi_eq_2/(3-2*self.Phi_eq_2) #para valor preciso considerar q!=0
                    #self.R_trans_mn = self.r_trans*np.sqrt(1-self.Phi_trans)
                    #self.fluid_trans.append(self.R_trans_mn)
                    #self.fluid_trans_mirror.append(-1*self.R_trans_mn)
                    #self.fluid_trans_H.append(self.L_out_2-2*self.trans+(self.R_out_2-self.r_trans)/(np.tan(self.fr.Trans_angle_2*np.pi/180)))
                
                ##### ######
                
                
                fig, ax4 = plt.subplots()
                ax4.grid()
                
                ax4.plot(self.inj_D_in,self.inj_H, linewidth = 2, color = 'red')
                ax4.plot(self.inj_D_in_mirror,self.inj_H, linewidth = 2, color = 'red')
                ax4.plot(self.inj_D_out,self.inj_H_out, linewidth = 2, color = 'red')
                ax4.plot(self.inj_D_out_mirror,self.inj_H_out, linewidth = 2, color = 'red')
                ax4.plot(self.inj_F,self.inj_F_H, linewidth = 2, color = 'red')
                ax4.plot(self.inj_F_mirror,self.inj_F_H, linewidth = 2, color = 'red')
                
                ax4.plot(self.inj_D_top,self.inj_H_top, linewidth = 2, color = 'red')
                
                ax4.plot(self.inj_1_Sfluid_D,self.inj_1_S_H, linewidth = 1, color = 'red', linestyle = '--')
                ax4.plot(self.inj_1_Sfluid_D_mirror,self.inj_1_S_H, linewidth = 1, color = 'red', linestyle = '--')
                ax4.plot(self.inj_1_Sfluid_Dtrans,self.inj_1_Sfluid_Htrans, linewidth = 1, color = 'red', linestyle = '--')
                ax4.plot(self.inj_1_Sfluid_Dtrans_mirror,self.inj_1_Sfluid_Htrans, linewidth = 1, color = 'red', linestyle = '--')
                
                
                ax4.plot(self.inj_1_S_D,self.inj_1_S_H, linewidth = 2, color = 'red')
                ax4.plot(self.inj_1_S_D_mirror,self.inj_1_S_H, linewidth = 2, color = 'red')
                ax4.plot(self.inj_1_Strans_D,self.inj_1_Strans_H, linewidth = 2, color = 'red')
                ax4.plot(self.inj_1_Strans_D_mirror,self.inj_1_Strans_H, linewidth = 2, color = 'red')
                
                ax4.plot(self.circ_x_1,self.inj_1_in_top, linewidth = 1, color = 'red')
                ax4.plot(self.circ_x_1,self.inj_1_in_bottom, linewidth = 1, color = 'red')
                
    
                ax4.plot(self.fluid_D_1,self.inj_H, linewidth = 1, color = 'red', linestyle = '--')
                ax4.plot(self.fluid_D_1_mirror,self.inj_H, linewidth = 1, color = 'red', linestyle = '--')
                ax4.plot(self.x_1,self.y_1, linewidth = 1, color = 'red', linestyle = '--')
                ax4.plot(self.x_1_mirror,self.y_1, linewidth = 1, color = 'red', linestyle = '--')
                
                ax4.plot(self.inj_D_2_top,self.inj_H_2_top, linewidth = 2, color = '0.3')
                ax4.plot(self.inj_D_2_top_mirror,self.inj_H_2_top, linewidth = 2, color = '0.3')
                
                
                ax4.plot(self.inj_2_Sfluid_D,self.inj_2_S_H, linewidth = 1, color = '0.3', linestyle = '--')
                ax4.plot(self.inj_2_Sfluid_D_mirror,self.inj_2_S_H, linewidth = 1, color = '0.3', linestyle = '--')
                ax4.plot(self.inj_2_Sfluid_Dtrans,self.inj_2_Sfluid_Htrans, linewidth = 1, color = '0.3', linestyle = '--')
                ax4.plot(self.inj_2_Sfluid_Dtrans_mirror,self.inj_2_Sfluid_Htrans, linewidth = 1, color = '0.3', linestyle = '--')
                
                
                ax4.plot(self.inj_2_S_D,self.inj_2_S_H, linewidth = 2, color = '0.3')
                ax4.plot(self.inj_2_S_D_mirror,self.inj_2_S_H, linewidth = 2, color = '0.3')
                ax4.plot(self.inj_2_Strans_D,self.inj_2_Strans_H, linewidth = 2, color = '0.3')
                ax4.plot(self.inj_2_Strans_D_mirror,self.inj_2_Strans_H, linewidth = 2, color = '0.3')
                
                ax4.plot(self.circ_x_2,self.inj_2_in_top, linewidth = 1, color = '0.3')
                ax4.plot(self.circ_x_2,self.inj_2_in_bottom, linewidth = 1, color = '0.3')
                
                ax4.plot(self.inj_2_D_in,self.inj_2_H, linewidth = 2, color = '0.3')
                ax4.plot(self.inj_2_D_in_mirror,self.inj_2_H, linewidth = 2, color = '0.3')
                #ax4.plot(self.inj_2_D_out,[0,self.L_n_2], linewidth = 2, color = '0.3')
                #ax4.plot(self.inj_2_D_out_mirror,[0,self.L_n_2], linewidth = 2, color = '0.3')
                ax4.plot(self.inj_2_F,self.inj_2_F_H, linewidth = 2, color = '0.3')
                ax4.plot(self.inj_2_F_mirror,self.inj_2_F_H, linewidth = 2, color = '0.3')
                
                ax4.plot(self.fluid_D_2,self.inj_2_H, linewidth = 1, color = '0.3', linestyle = '--')
                ax4.plot(self.fluid_D_2_mirror,self.inj_2_H, linewidth = 1, color = '0.3', linestyle = '--')
                #ax4.plot(x_2,y_2, linewidth = 1, color = 'blue', linestyle = '--')
                #ax4.plot(x_2_mirror,y_2, linewidth = 1, color = 'blue', linestyle = '--')
                
                
                ax4.plot(self.inj_3_D_in,self.inj_3_H, linewidth = 2, color = '0.3')
                ax4.plot(self.inj_3_D_in_mirror,self.inj_3_H, linewidth = 2, color = '0.3')
                ax4.plot(self.inj_3_Dtrans_in,self.inj_3_Htrans, linewidth = 2, color = '0.3')
                ax4.plot(self.inj_3_Dtrans_in_mirror,self.inj_3_Htrans, linewidth = 2, color = '0.3')
                
                ax4.plot(self.fluid_D_3,self.inj_3_H, linewidth = 1, color = '0.3', linestyle = '--')
                ax4.plot(self.fluid_D_3_mirror,self.inj_3_H, linewidth = 1, color = '0.3', linestyle = '--')
                ax4.plot(self.fluid_trans,self.fluid_trans_H, linewidth = 1, color = '0.3', linestyle = '--')
                ax4.plot(self.fluid_trans_mirror,self.fluid_trans_H, linewidth = 1, color = '0.3', linestyle = '--')
                
                ax4.plot(self.x_3,self.y_3, linewidth = 1, color = '0.3', linestyle = '--')
                ax4.plot(self.x_3_mirror,self.y_3, linewidth = 1, color = '0.3', linestyle = '--')
                
                
                #ax4.plot(self.inj_23_D_in,self.inj_23_H, linewidth = 2, color = '0.3')
                #ax4.plot(self.inj_23_D_in_mirror,self.inj_23_H, linewidth = 2, color = '0.3')
                
                ax4.set_xlabel("[mm]")
                ax4.set_ylabel("[mm]")
                
                self.min_x = min([-1*(self.R_out_2+self.fr.t_w+2.0),-1*(self.R_s_1+2.0)])
                self.max_x = max([(self.R_out_2+self.fr.t_w+2.0),(self.R_s_1+2.0)])
                self.min_y = -5.0
                self.max_y = self.H_1_s_top+5.0
                
                #ax4.axis(xlim=(self.min_x, self.max_x), ylim=(self.min_y, self.max_y))
                ax4.axis('equal')
                
        #####  end of Bipropellant only  #####

        if(self.fr.Config == "Mono\n"):
            
            ##### ##### #####
            
            self.H_1_trans = (self.R_s_1-self.R_n_1)/np.tan(self.fr.Trans_angle_1*np.pi/180)
            self.inj_D_top = [-1*self.R_s_1,self.R_s_1]
            self.inj_H_top = [self.L_n_1 + self.H_1_trans + self.L_s_1,self.L_n_1 + self.H_1_trans + self.L_s_1]
            self.fluid_D_1 = [self.r_mn_1,self.r_mn_1]
            self.fluid_D_1_mirror = [-1*self.r_mn_1,-1*self.r_mn_1]
            
            self.inj_H = [0.0,self.L_n_1]
            self.inj_D_in = [self.R_n_1,self.R_n_1]
            self.inj_D_in_mirror = [-1*self.R_n_1,-1*self.R_n_1]
            
            self.inj_1_Strans_D = [self.R_n_1,self.R_s_1]
            self.inj_1_Strans_D_mirror = [-1*self.R_n_1,-1*self.R_s_1]
            self.inj_1_Strans_H = [self.L_n_1,self.L_n_1 + self.H_1_trans]
            
            self.inj_1_Sfluid_D = [self.r_mk_1,self.r_mk_1]
            self.inj_1_Sfluid_D_mirror = [-1*self.r_mk_1,-1*self.r_mk_1]
            
            self.inj_1_S_D = [self.R_s_1,self.R_s_1]
            self.inj_1_S_D_mirror = [-1*self.R_s_1,-1*self.R_s_1]
            self.inj_1_S_H = [self.L_n_1 + self.H_1_trans,self.L_n_1 + self.H_1_trans + self.L_s_1]
            
            self.inj_1_Sfluid_Dtrans = [self.r_mn_1,self.r_mk_1]
            self.inj_1_Sfluid_Dtrans_mirror = [-1*self.r_mn_1,-1*self.r_mk_1]
            self.inj_1_Sfluid_Htrans = [self.L_n_1,self.L_n_1 + self.H_1_trans]
            
            self.inj_1_in_top = []
            self.inj_1_in_bottom = []
            
            self.circ_x_1 = np.arange(-1*self.R_s_1,-1*self.R_s_1+2*self.r_in_1+0.005,0.005)
            
            for x in self.circ_x_1:
                self.inj_1_in_top.append((self.L_n_1 + self.H_1_trans + self.L_s_1-self.r_in_1) + np.sqrt(self.r_in_1**2-(x+self.R_s_1-self.r_in_1)**2))
                self.inj_1_in_bottom.append((self.L_n_1 + self.H_1_trans + self.L_s_1-self.r_in_1) - np.sqrt(self.r_in_1**2-(x+self.R_s_1-self.r_in_1)**2))
            
            self.y_i_1 = 0
            self.y_f_1 = -5
            self.y_step_1 = (self.y_f_1-self.y_i_1)/100
            self.y_1 = np.arange(self.y_i_1,self.y_f_1+self.y_step_1,self.y_step_1)
            
            self.x_1 = []
            self.x_1_mirror = []
            
            self.R_av_1 = (self.R_n_1 + self.r_mn_1)/2
            
            for y in self.y_1:
                self.x_temp = self.R_av_1*np.sqrt(1+(y-self.fr.recess)**2/(self.R_av_1*np.tan((90-0.5*self.Alpha2_eq_1_calc)*np.pi/180))**2)
                self.x_1.append(self.x_temp)
                self.x_1_mirror.append(-1*self.x_temp)   
             
            
            
            ##### ##### #####
            
            fig, ax5 = plt.subplots()
            ax5.grid()
        
            ax5.plot(self.inj_D_in,self.inj_H, linewidth = 2, color = 'red')
            ax5.plot(self.inj_D_in_mirror,self.inj_H, linewidth = 2, color = 'red')            
            ax5.plot(self.inj_D_top,self.inj_H_top, linewidth = 2, color = 'red')
            
            ax5.plot(self.inj_1_Sfluid_D,self.inj_1_S_H, linewidth = 1, color = 'red', linestyle = '--')
            ax5.plot(self.inj_1_Sfluid_D_mirror,self.inj_1_S_H, linewidth = 1, color = 'red', linestyle = '--')
            ax5.plot(self.inj_1_Sfluid_Dtrans,self.inj_1_Sfluid_Htrans, linewidth = 1, color = 'red', linestyle = '--')
            ax5.plot(self.inj_1_Sfluid_Dtrans_mirror,self.inj_1_Sfluid_Htrans, linewidth = 1, color = 'red', linestyle = '--')
            
            
            ax5.plot(self.inj_1_S_D,self.inj_1_S_H, linewidth = 2, color = 'red')
            ax5.plot(self.inj_1_S_D_mirror,self.inj_1_S_H, linewidth = 2, color = 'red')
            ax5.plot(self.inj_1_Strans_D,self.inj_1_Strans_H, linewidth = 2, color = 'red')
            ax5.plot(self.inj_1_Strans_D_mirror,self.inj_1_Strans_H, linewidth = 2, color = 'red')
            
            ax5.plot(self.circ_x_1,self.inj_1_in_top, linewidth = 1, color = 'red')
            ax5.plot(self.circ_x_1,self.inj_1_in_bottom, linewidth = 1, color = 'red')
            

            ax5.plot(self.fluid_D_1,self.inj_H, linewidth = 1, color = 'red', linestyle = '--')
            ax5.plot(self.fluid_D_1_mirror,self.inj_H, linewidth = 1, color = 'red', linestyle = '--')
            ax5.plot(self.x_1,self.y_1, linewidth = 1, color = 'red', linestyle = '--')
            ax5.plot(self.x_1_mirror,self.y_1, linewidth = 1, color = 'red', linestyle = '--')
            
            
            ax5.set_xlabel("[mm]")
            ax5.set_ylabel("[mm]")
            
            #self.min_x = min([-1*(self.R_out_2+self.fr.t_w+2.0),-1*(self.R_s_1+2.0)])
            #self.max_x = max([(self.R_out_2+self.fr.t_w+2.0),(self.R_s_1+2.0)])
            #self.min_y = -5.0
            #self.max_y = self.H_1_s_top+5.0
            
            #ax5.axis(xlim=(self.min_x, self.max_x), ylim=(self.min_y, self.max_y))
            ax5.axis('equal')
            