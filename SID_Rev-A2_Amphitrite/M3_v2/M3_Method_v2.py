
"""
Created on Mon Sep  7 07:41:29 2020

@author: A. Goulart
"""

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt

from M3_FileReader_v2 import FileReader_3


class Method_3:
    def __init__(self, foldername, filename):
        self.it_lim = 1000
        self.erro_max = 1e-4
        self.sample_size = 1000
        
        self.delta_p_i = 0.05 #bar
        self.delta_p_f = 10.0 #bar
        self.delta_p_step = 0.01
        
        
        self.fr = FileReader_3()
        self.fr.setFolderName(foldername)
        self.fr.setFileName(filename)
        
        self.fmt = '| {{:^{}s}} | {{:^{}s}} |'.format(25, 15) # widths only

    
    def set_delta_p_i(self,n):
        if(n == 0):
            self.delta_p_i = 0.01
        else:
            self.delta_p_i = n
        
    def get_delta_p_i(self):
        print("delta_p_i = ", self.delta_p_i)


    def set_delta_p_f(self,n):
        self.delta_p_f = n
        
    def get_delta_p_f(self):
        print("delta_p_f = ", self.delta_p_f)


    def set_delta_p_step(self,n):
        self.delta_p_step = n
        
    def get_delta_p_step(self):
        print("delta_p_step = ", self.delta_p_step)


    def set_it_lim(self,n):
        self.it_lim = n
        
    def get_it_lim(self):
        print("it_lim = ", self.it_lim)

        
    def set_erro_max(self,n):
        self.erro_max = n
        
    def get_erro_max(self):
        print("erro_max = ", self.erro_max)
  
    
    def set_sample_size(self,n):
        self.sample_size = n
        
    def get_sample_size(self):
        print("sample_size = ", self.sample_size)
  
    
    def print_data(self, label, value):
        print(self.fmt.format(str(label), str(value)))



    ##### Métodos para cálculo do injetore #####
    
    def calc_phi(self, geom_char):
        self.coefs = [geom_char, -2, 4, -2]
        self.roots = np.roots(self.coefs)

        self.count = 0
        
        self.phi_eq = -1
        
        for root in self.roots:
            if(0 < root < 1 and root.imag == 0):
                self.count+=1
                self.phi = root.real
        
        if(self.phi == -1):
            print("Error in calc_phi:\n No valid root found.")
            
        if(self.count > 1):
            print("Error in calc_phi:\n Multiple valid values found")
            print(self.roots)
        
        return self.phi
    
    
    
    def calc_geom_char(self,r_in_pos, r_n, n, r_in_orf, friction_coef):
        return r_in_pos*1e-3*r_n*1e-3/(n*(r_in_orf*1e-3)**2+friction_coef/2*r_in_pos*1e-3*(r_in_pos*1e-3-r_n*1e-3))
    
    def calc_2alpha(self, phi):
        return 2*180/np.pi*np.arctan(np.sqrt(2*(1-phi)/phi))
    
    def calc_reynolds_in(self,m_in,din_visc,n,r_in_orf):
        return (2*m_in)/(din_visc*n*np.pi*r_in_orf)
    
    def calc_friction_coef(self, reynolds_in):
        #return 10**(25.8/(np.log(reynolds_in))-2)  #Bayvel
        return 0.3164*reynolds_in**(-0.25)          #Bazarov
    
    
    ## Calcula o valor da perda hidráulica no injetor (Ksi)
    def calc_hyd_loss(self, inlet_type, r_s, l_in, inlet_radius, r_in_orf, friction_coef):
        if(inlet_type == "curved\n"):
            self.alpha_in = 90 - 180/np.pi * np.arctan(r_s/l_in)
            self.ksi_in = -1*self.alpha_in/150 + 1.1
            
        elif(inlet_type == "straigth\n"):
            self.ksi_in = 0.5*np.exp(-1.4*(inlet_radius/(2*r_in_orf)))
            
        else:
            print("Error in calc_hyd_loss:\n Inlet type not recognized")
        
        return self.ksi_in + friction_coef*l_in/(2*r_in_orf)
    
    ## calcula o parâmetro geométrico característico corrigido (A_eq) [Bazarov, eq(100)]
    def calc_geom_char_eq(self, r_in_pos, r_n, n, r_in_orf, friction_coef):
        return r_in_pos*r_n/(n*r_in_orf**2+friction_coef/2*r_in_pos*(r_in_pos-r_n))
    
    
    ## Calcula o coeficiente de descarga equivalente [Bazarov, eq(99)]
    def calc_mu_eq(self, phi_eq):
        return phi_eq**1.5/np.sqrt(2-phi_eq)
    
    ## Calcula o coeficiente de descarga real do injetor [Bazarov, eq(99)]
    def calc_mu_i(self, mu_eq, geom_char, r_in_pos, r_n, ksi):
        return mu_eq/np.sqrt(1+ksi*mu_eq*(geom_char*r_in_pos/r_n)**2)

    
    ##### Métodos - Cálculos adicionais #####
    
    ##Calcula o valor da não-uniformidade percentual esperada do injetor [Bayvel]
    def calc_non_uniformity(self, r_in_pos, r_n, n, l_s, r_s, phi_eq):
        return 23.7/((r_in_pos/r_n)**2.7*n**1.34*phi_eq**1.1*(l_s/(2*r_s))**0.15)
    
    
    
    ##### Cálculo do injetor #####
    
    def run_M3(self):
            
        self.fr.read()
        
        
        #SIG - Swirl Injector from Geometry
        print("\n\n Method M3 \n")

        ##### ##### ##### ##### ##### ##### #####
        ##### #####  ST1 Calculation  ##### #####
        ##### ##### ##### ##### ##### ##### #####
        
        
        print("m_target: \t%.2f \tg/s\n"%(np.round(self.fr.m_0,2)))

        print("\n ST1 - Geometria do injetor:")
 
        self.print_data("n", "%.1f"%(self.fr.n))

        self.print_data("r_in_orf", "%.2f mm"%(self.fr.r_in_orf))

        self.print_data("l_in", "%.2f mm"%(self.fr.l_in))

        self.print_data("r_s", "%.2f mm"%(self.fr.r_s))

        self.print_data("l_s", "%.2f mm"%(self.fr.l_s))
        
        self.print_data("r_n", "%.2f mm"%(self.fr.r_n))
        
        self.print_data("r_in_pos", "%.2f mm"%(self.fr.r_in_pos))
        
        
        self.n_conv = []
        self.m_ideal = []
        self.m_real = []

        self.pressure_range = np.arange(self.delta_p_i, self.delta_p_f+self.delta_p_step, self.delta_p_step)

        for self.delta_p in self.pressure_range:
            
            print("\n\nCalculating for delta_p = ", self.delta_p, " bar")
            
            self.erro_m = self.erro_max + 1
            self.m_prev = self.erro_max + 1
            self.contador = 0
            
            while self.erro_m >= self.erro_max and self.contador < self.it_lim:
                
                ##### ##### ##### ##### ##### ###### #####
                ##### ST1 - Cálculo do Estágio Ideal #####
                ##### ##### ##### ##### ##### ###### #####
                
                
                ### Ideal Injector Parameters  ###
                
                self.A = self.calc_geom_char(self.fr.r_in_pos, self.fr.r_n, self.fr.n, self.fr.r_in_orf, 0) 
                 
                
                self.phi = self.calc_phi(self.A)
                

                #Coeficiente de Descarga ideal
                self.mu = self.phi*np.sqrt(self.phi/(2-self.phi))
                
                #raio da camada líquida ideal
                self.r_mn = self.fr.r_n*np.sqrt(1-self.phi)
                
                #espessura da camada líquida ideal
                self.t_mn = self.fr.r_n-self.r_mn
                

                #Vazão mássica ideal
                self.m = self.mu*(np.pi*self.fr.r_n**2)*1e-6*np.sqrt(2*self.delta_p*1e5*self.fr.rho)*1e3 #vazão mássica ideal [g/s]
                
                #Coeficiente de redução de área do estágio
                self.C = self.fr.r_in_pos/(self.fr.r_n)
                
                ##### Cálulo do ângulo de spray ideal #####
                
                self.alpha2 = self.calc_2alpha(self.phi)
        
        
        
                
                ##### ##### ##### ##### ##### ##### #####
                ##### ST1 - Cálculo do Estágio Real #####
                ##### ##### ##### ##### ##### ##### #####
                
                
                #Número de Reynolds nos canais tangenciais
                if(self.contador == 0):
                    self.rey_in = self.calc_reynolds_in(self.m,self.fr.din_visc,self.fr.n,self.fr.r_in_orf)
                else:
                    self.rey_in = self.calc_reynolds_in(self.m_prev,self.fr.din_visc,self.fr.n,self.fr.r_in_orf)
                print("\nrey_in = ", self.rey_in)
                
                #Coeficiente de resistência de superfície
                self.friction_coef = self.calc_friction_coef(self.rey_in)
                        
                
                #Parâmetro Geométrico equivalente
                self.A_eq = self.calc_geom_char(self.fr.r_in_pos, self.fr.r_n, self.fr.n, self.fr.r_in_orf, self.friction_coef) 
                        
                
                ##### Cálculo do Coeficiente de área livre real
                self.phi_eq = self.calc_phi(self.A_eq)
                
        
                #raio da camada líquida real
                self.r_mn_eq = self.fr.r_n*np.sqrt(1-self.phi_eq)
        
                #Espessura de camada líquida real
                self.t_mn_eq = self.fr.r_n*(1-np.sqrt(1-self.phi_eq)) #mm
                
                
                
                #Coeficiente de perda de momento angular
                self.K = self.A_eq/self.A
                
                
                
                ##### Cálculo da perda nos canais tangenciais #####
        
                self.ksi = self.calc_hyd_loss(self.fr.inlet_type, self.fr.r_s, self.fr.l_in, self.fr.inlet_radius, self.fr.r_in_orf, self.friction_coef)
                
                
                #Coeficiente de Descarga
                self.mu_eq = self.calc_mu_eq(self.phi_eq)

                self.mu_i = self.calc_mu_i(self.mu_eq, self.A_eq, self.fr.r_in_pos, self.fr.r_n, self.ksi)
                
                #Vazão mássica real
                self.m_i = np.pi*(self.fr.r_n*1e-3)**2*np.sqrt(2*self.fr.rho*self.delta_p*1e5)/np.sqrt((2-self.phi_eq)/self.phi_eq**3+self.ksi*self.A**2/(self.fr.r_in_pos)**2)*1e3
                print("m_i = ", self.m_i, " g/s")
                
                self.erro_m = abs((self.m_i-self.m_prev)/self.m_prev)
                print("erro_m = ", self.erro_m)
                
                self.m_prev = self.m_i
                
                ##### Cálulo do ângulo de spray real #####
                self.alpha2_eq = self.calc_2alpha(self.phi_eq)
                    
                ##### Cálculo da não-uniformidade #####
                self.I = self.calc_non_uniformity(self.fr.r_in_pos, self.fr.r_n, self.fr.n, self.fr.l_s, self.fr.r_s, self.phi_eq)
                
                
                self.contador +=1
                print("contador = ", self.contador)
                
            self.n_conv.append(self.contador)
            
            self.m_ideal.append(self.m)
            self.m_real.append(self.m_i)
            
            ##### Tests #####
            if(self.t_mn <= 0): print("\nERROR: invalid t_mn value (delta_p = %.2f bar)"%(self.delta_p))
            elif(self.t_mn <= 0.1): print("\nCaution! Low value of t_mn (delta_p = %.2f bar)"%(self.delta_p))
                
            if(self.t_mn_eq <= 0): print("\nERROR: invalid t_mn_eq value (delta_p = %.2f bar)"%(self.delta_p))
            elif(self.t_mn_eq <= 0.1): print("\nCaution! Low value of t_mn_eq (delta_p = %.2f bar)"%(self.delta_p))
            
            if(self.rey_in < 2600): print("\nERROR: Low Reynolds number (%.2E < 2600) (delta_p = %.2f bar)"%(self.rey_in, self.delta_p))
            elif(self.rey_in > 5e4): print("\nCaution: High Reynolds number (%.2E > 5e4) (delta_p = %.2f bar)"%(self.rey_in, self.delta_p))
                
            
                        
                
            if(self.contador >= self.it_lim):
                print ("\n\tWarning:\n iteration limit reached (delta_p = %.2f bar)"%(self.delta_p))
                print ("\terror = ", self.erro_m)
            
            
            
        self.step_x = np.ceil(self.contador/10)
        
        fig, (ax1) = plt.subplots()
        fig.suptitle("Solution convergence")
        
        ax1.plot(self.pressure_range,self.n_conv)
        ax1.set_xlabel("delta_p [bar]")
        ax1.set_ylabel("# interations")
        #ax1.set_xticks(np.arange(0,self.contador+self.step_x,self.step_x))
        ax1.grid()
        
        fig, ax11 = plt.subplots(figsize=(14, 7))
        ax11.grid()
        ax11.set_title("ṁ x ΔP") 
        ax11.plot(self.pressure_range, self.m_ideal, 'k--', markersize=1, label = "ṁ fluido ideal", color = '0.5')
        ax11.plot(self.pressure_range, self.m_real, 'k:', markersize=3, label = "ṁ fluido viscoso")
        ax11.plot(self.pressure_range, [self.fr.m_0]*len(self.pressure_range), 'k-.', markersize=1, label = "ṁ alvo")
        ax11.plot(self.fr.delta_p0, self.fr.m_0, '.', markersize=20, label = "ṁxΔP alvo", color = '0.5')
        ax11.legend()
        ax11.set_ylabel("ṁ [g/s]")
        ax11.set_xlabel("ΔP [bar]")
        #ax11.set_xticks(np.arange(0.0,10.0+0.2,0.2))
        #self.y_min = round(min([min(self.Mass_ideal),min(self.Mass_real),min(self.M_target),self.fr.m_1]),0)
        #self.y_max = round(max([max(self.Mass_ideal),max(self.Mass_real),max(self.M_target),self.fr.m_1]),0)
        #self.y_gap = (self.y_max - self.y_min)/20
        #ax11.set_yticks(np.arange(self.y_min,self.y_max+self.y_gap,self.y_gap))
        ax11.tick_params(axis='x', labelsize=8)
        ax11.tick_params(axis='y', labelsize=8)
        #ax11.set_xlim(0.0,10.0)
        #ax11.grid()
        
        
        
        
        
        
        
        
   
    