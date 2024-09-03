# -*- coding: utf-8 -*-
"""
@author: A. Goulart
"""

import numpy as np

import matplotlib.pyplot as plt

from M1_FileReader import FileReader_1



class Method_1:
    def __init__(self, foldername, filename):
        self.it_lim = 1000
        self.erro_max = 1e-4
        
        self.fontsize = 10
        self.offset = 26
        self.dpi_value = 800
        self.show_it = 1
        self.fmt = '| {{:^{}s}} | {{:^{}s}} |'.format(25, 15) # widths only
        
        self.fr = FileReader_1()
        self.fr.setFolderName(foldername)
        self.fr.setFileName(filename)

        plt.isinteractive()


    def set_show_it(self,show_val):
        self.show_it = show_val
        
    def get_show_it(self):
        return self.show_it
        
    def set_it_lim(self,n):
        self.it_lim = n
        
    def get_it_lim(self):
        return self.it_lim
        
    def set_erro_max(self,n):
        self.erro_max = n
        
    def get_erro_max(self):
        return self.erro_max
        
    def set_textbox_fontsize(self,fs):
        self.fontsize = int(fs)
    
    def get_textbox_fontsize(self):
        return self.fontsize
        
    def set_textbox_offset(self,off):
        self.offset = off
        
    def get_textbox_offset(self):
        return self.offset
        
    def show_plot_inj(self,answer):
        self.sh_plot_inj = answer
        
    def show_plot_error(self,answer):
        self.sh_plot_error = answer
        
    def set_dpi(self,n):
        self.dpi_value = n
        
    def get_dpi(self):
        return self.dpi_value
        
    def print_data(self, label, value):
        print(self.fmt.format(str(label), str(value)))
    
    def check_value(self, var_name, min_value, max_value, var_value):
        if(var_value > max_value):
            print("ERROR: %s value out of bounds [%.2f, %.2f] (current value: %.2f)"(var_name,min_value,max_value,var_value))
            
        elif(var_value < min_value):
            print("ERROR: %s value out of bounds [%.2f, %.2f] (current value: %.2f)"(var_name,min_value,max_value,var_value))
            
    
    ##### Métodos - Injetor Ideal #####
    
    def calc_phi(self, alpha2):
        return 2/(2+(np.tan(alpha2/2*np.pi/180))**2)
        '''
        self.gamma_phi = np.tan(alpha2/2*np.pi/180)
        
        self.a_phi = 2*(1 + self.gamma_phi**2)
        self.b_phi = self.gamma_phi**2
        self.c_phi = -self.gamma_phi**2
        
        self.delta_phi = self.b_phi**2 - 4*self.a_phi*self.c_phi#(9*self.gamma_phi + 1) / (4*self.gamma_phi)
        
        if(self.delta_phi < 0):
            print("\nError in Calc_Phi_eq:\n No valid root found.\n")
        
        elif(self.delta_phi > 0):
            print("\nWarning! Multiple valid values found in Calc_Phi")
            self.calc_phi_1 = 1- (-self.b_phi + np.sqrt(self.delta_phi)/(2*self.a_phi))
            self.calc_phi_2 = 1- (-self.b_phi - np.sqrt(self.delta_phi)/(2*self.a_phi))
            print("phi_1 = %.2f \nphi_2 = %.2f\n"%(self.calc_phi_1, self.calc_phi_2))
            return self.calc_phi_1
            
        else:
            return 1 - (-self.b_phi / (2*self.a_phi))
        '''
    
    def calc_geom_char(self, phi):
        return (1-phi)*np.sqrt(2/phi)/phi
        
    def calc_mu(self, phi):
        return phi*np.sqrt(phi/(2-phi))
    
    def calc_a(self,phi):
        return 2*((1-phi)**2)/(2-phi)
    
    def calc_2alpha(self, coef, coef_type):
        if(coef_type == "phi"):
            return 2*180/np.pi*np.arctan(np.sqrt(2*(1-coef)/coef))
        elif(coef_type == "a"):
            return 2*180/np.pi*np.arctan(coef/(1-coef))
        else:
            print("ERROR: input not understood in calc_2alpha")
            
    def calc_2alpha_avg(self, mu, A, phi, ksi, C):
        return 2*(180/np.pi)*np.arcsin(2*mu*A/((1+np.sqrt((1-phi))*(1-ksi*mu**2*(A/(C))**2))))
    
    ##### Métodos - Injetor Real #####
    
    ## Calcula o valor do nº de Reynolds nos canais tangenciais do injetor
    def calc_reynolds_in(self,m_in,din_visc,n,r_in_orf,rho):
        return (2*m_in*1e-3)/(din_visc*n*np.pi*(r_in_orf*1e-3))
        #return (0.637*m_in)/(np.sqrt(n)*r_in_orf*rho*din_visc)
    
    ## calcula o valor do coeficiente de atrito (Lambda)
    def calc_friction_coef(self, reynolds_in):
        #return 10**(25.8/(np.log(reynolds_in))-2)  #Bayvel
        return 0.3164*reynolds_in**(-0.25)          #Bazarov
    
    ## Calcula o valor da perda hidráulica no injetor (Ksi)
    def calc_hyd_loss(self, inlet_type, r_s, l_in, inlet_radius, r_in_orf, friction_coef):
        if(inlet_type == "curved\n"):
            self.alpha_in = 90 - 180/np.pi * np.arctan(r_s/l_in) #eq(105)
            self.ksi_in = -1*self.alpha_in/150 + 1.1
            
        elif(inlet_type == "straight\n"):
            self.ksi_in = 0.5*np.exp(-1.4*(inlet_radius/(2*r_in_orf)))
            
        else:
            print("Error in calc_hyd_loss:\n Inlet type not recognized")
        
        return self.ksi_in + friction_coef*l_in/(2*r_in_orf)

    ## calcula o parâmetro geométrico característico corrigido (A_eq) [Bazarov, eq(100)]
    def calc_geom_char_eq(self, r_in_pos, r_n, n, r_in_orf, friction_coef):
        return (r_in_pos*r_n)/(n*r_in_orf**2+(friction_coef/2)*r_in_pos*(r_in_pos-r_n))
    
   
    ## Calcula o Valor do Coeficiente de Preenchimento corrigido (Phi_eq) [Bazarov, eq(99), eq(100)]
    def calc_phi_eq(self, geom_char_eq):
        self.coefs = [geom_char_eq**2, -2, 4, -2]
        self.roots = np.roots(self.coefs)

        self.count = 0
        
        self.phi_eq = -1
        
        for root in self.roots:
            if(0 < root < 1 and root.imag == 0):
                self.count+=1
                self.phi_eq = root.real
        
        if(self.phi_eq == -1):
            print("Error in Calc_Phi_eq:\n No valid root found.")
            
        if(self.count > 1):
            print("Error in Calc_Phi_eq:\n Multiple valid values found")
            print(self.roots)
        
        return self.phi_eq
    
    ## Calcula o coeficiente de descarga equivalente [Bazarov, eq(99)]
    def calc_mu_eq(self, phi_eq):
        return phi_eq**1.5/np.sqrt(2-phi_eq)
    
    ## Calcula o coeficiente de descarga real do injetor [Bazarov, eq(99)]
    def calc_mu_i(self, mu_eq, geom_char, r_in_pos, r_n, ksi):
        return mu_eq/np.sqrt(1+ksi*mu_eq**2*(geom_char/(r_in_pos/r_n))**2)
    
    ##### Métodos - Cálculos adicionais #####
    
    ##Calcula o valor da não-uniformidade percentual esperada do injetor [Bayvel]
    def calc_non_uniformity(self, r_in_pos, r_n, n, l_s, r_s, phi_eq):
        return 23.7/((r_in_pos/r_n)**2.7*n**1.34*phi_eq**1.1*(l_s/(2*r_s))**0.15)
    
    
    ##### Roda o método de dimensionamento do injetor #####
    
    def run_M1(self):
        
        self.fr.read()
        
        print("Config: ", self.fr.config)

        ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
        ##### ##### ##### Parte 1 - Cálculo do ST1 #### ##### ##### ##### ##### #####
        ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
        
        self.phi_1 = self.calc_phi(self.fr.alpha2_1)
        
        self.A_1 = self.calc_geom_char(self.phi_1)
        
        self.mu_1 = self.calc_mu(self.phi_1)
        
        #self.r_n_1 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(self.fr.m_1*1e-3/(self.mu_1*np.sqrt(2*self.fr.rho_1*self.fr.delta_p_1*1e5)))
        #self.r_n_1 = 1e3 * np.sqrt(self.fr.m_1*1e-3/(np.pi*self.mu_1*np.sqrt(2*self.fr.rho_1*self.fr.delta_p_1*1e5)))
        self.r_n_1 = 1e3*0.475*np.sqrt(self.fr.m_1*1e-3 / (self.mu_1 * np.sqrt(self.fr.rho_1*self.fr.delta_p_1*1e5)))
        
        print("\n\n ST1 Initial Values:")
        self.print_data("A_1", "%.2f"%(self.A_1))
        self.print_data("Phi_1", "%.3f"%(self.phi_1))
        self.print_data("R_n_1", "%.3f mm"%(self.r_n_1))
        self.print_data("2Alpha", "%.1f deg"%(self.fr.alpha2_1))

        
        self.inlet_diameter = []
        self.outlet_error = []
        
        self.erro_r = self.erro_max + 1
        self.contador = 0
        
        if(self.show_it == 1): print("\n ST1 Iteration data:")
        while self.erro_r >= self.erro_max and self.contador < self.it_lim:
            
            
            ##### Iteração do injetor ideal #####
            
            self.r_in_pos_1 = self.fr.opening_coef_1*self.r_n_1
            
            self.r_in_orf_1 = np.sqrt(self.r_in_pos_1*self.r_n_1/(self.fr.n_1*self.A_1))
            
            self.inlet_diameter.append(2*self.r_in_orf_1)
            
            self.r_s_1 = self.r_in_pos_1+self.r_in_orf_1
            
            self.l_in_1 = self.fr.ratio_l_in_1*2*self.r_in_orf_1
            
            self.l_n_1 = self.fr.ratio_l_n_1*2*self.r_n_1

            self.l_s_1 = self.fr.ratio_l_s_1*2*self.r_s_1
            
            self.a_1 = self.calc_a(self.phi_1)
                       
            self.alpha2_1 = self.calc_2alpha(self.phi_1, "phi")
            
            self.alpha2_e_1 = self.calc_2alpha(self.a_1, "a")
            
            ##### Iteração do injetor real #####
            
            self.rey_in_1 = self.calc_reynolds_in(self.fr.m_1, self.fr.din_visc_1, self.fr.n_1, self.r_in_orf_1,self.fr.rho_1)
            
            self.lambda_1 = self.calc_friction_coef(self.rey_in_1)
            
            self.ksi_1 = self.calc_hyd_loss(self.fr.in_type_1, self.r_s_1, self.l_in_1, self.fr.inlet_radius_1, self.r_in_orf_1, self.lambda_1)
            self.check_value("ksi_1", 0, float("inf"), self.ksi_1)
            
            self.A_eq_1 = self.calc_geom_char_eq(self.r_in_pos_1, self.r_n_1, self.fr.n_1, self.r_in_orf_1, self.lambda_1)
            
            self.phi_eq_1 = self.calc_phi_eq(self.A_eq_1)
            
            self.mu_eq_1 = self.calc_mu_eq(self.phi_eq_1)
            
            self.mu_i_1 = self.calc_mu_i(self.mu_eq_1, self.A_eq_1, self.r_in_pos_1, self.r_n_1, self.ksi_1)
            
            self.a_eq_1 = self.calc_a(self.phi_eq_1)
            
            self.alpha2_eq_1 = self.calc_2alpha(self.phi_eq_1, "phi")
            
            self.alpha2_e_eq_1 = self.calc_2alpha(self.a_eq_1, "a")
            
            ## cálculo do novo raio do bocal
            
            self.r_n_old = self.r_n_1
            
            #self.r_n_1 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(self.fr.m_1*1e-3/(self.mu_i_1*np.sqrt(2*self.fr.rho_1*self.fr.delta_p_1*1e5)))
            #self.r_n_1 = 1e3 * np.sqrt(self.fr.m_1*1e-3/(np.pi*self.mu_i_1*np.sqrt(2*self.fr.rho_1*self.fr.delta_p_1*1e5)))
            self.r_n_1 = 1e3*0.475*np.sqrt(self.fr.m_1*1e-3 / (self.mu_i_1 * np.sqrt(self.fr.rho_1*self.fr.delta_p_1*1e5)))
            
            self.erro_r = abs(self.r_n_old-self.r_n_1)/self.r_n_old
            
            self.outlet_error.append(self.erro_r)
            
            
            self.A_1 = self.r_in_pos_1*self.r_n_1/(self.fr.n_1*self.r_in_orf_1**2)
            
            self.phi_1 = self.calc_phi_eq(self.A_1)
            
            self.K_1 = self.A_eq_1/self.A_1
            self.check_value("K_1", 0, 1, self.K_1)
            
            self.r_mn_1 = self.r_n_1*np.sqrt(1-self.phi_eq_1)
            
            self.r_mk_1 = self.r_mn_1 * np.sqrt(2*(1-self.phi_eq_1)/(2-self.phi_eq_1))
            
            self.t_fluid_1 = self.r_n_1-self.r_mn_1
            
            self.C_1 = self.r_in_pos_1/self.r_n_1
            
            self.alpha2_avg_1 = self.calc_2alpha_avg(self.mu_i_1, self.A_eq_1, self.phi_eq_1, self.ksi_1, self.C_1)
            
            
            ## print dos dados da iteração
            if(self.show_it == 1):
                print("\nIt. %d)"%(self.contador+1))
                self.print_data("A_1", "%.4f"%(self.A_1))
                self.print_data("A_eq_1", "%.4f"%(self.A_eq_1))
                self.print_data("K_1", "%.4f"%(self.K_1))
                self.print_data("Phi_1", "%.4f"%(self.phi_1))
                self.print_data("Phi_eq_1", "%.4f"%(self.phi_eq_1))
                self.print_data("a_1", "%.4f"%(self.a_1))
                self.print_data("a_eq_1", "%.4f"%(self.a_eq_1))
                self.print_data("rey_in_1", "%.4f"%(self.rey_in_1))
                self.print_data("ksi_1", "%.4f"%(self.ksi_1))
                self.print_data("lambda_1", "%.4f"%(self.lambda_1))
                self.print_data("mu_i", "%.4f"%(self.mu_i_1))
                self.print_data("C", "%.4f"%(self.C_1))
                #self.print_data("2Alpha ideal", "%.2f deg"%(self.alpha2_1))
                #self.print_data("2Alpha real", "%.2f deg"%(self.alpha2_eq_1))
                #self.print_data("2Alpha ideal jusante", "%.2f deg"%(self.alpha2_e_1))
                #self.print_data("2Alpha real jusante", "%.2f deg"%(self.alpha2_e_eq_1))
                self.print_data("2Alpha_avg","%.2f deg"%(self.alpha2_avg_1))
                self.print_data("R_n_1", "%.3f mm"%(self.r_n_1))
                self.print_data("r_in_1", "%.3f mm"%(self.r_in_orf_1))
                self.print_data("Error", "%.2E"%(self.erro_r))
                
            self.contador += 1
        
        print("\n\n")
        
        ##Cálculo da não-uniformidade esperada do injetor
        self.I_1 = self.calc_non_uniformity(self.r_in_pos_1, self.r_n_1, self.fr.n_1, self.l_s_1, self.r_s_1, self.phi_eq_1)
              
        
        if(self.contador >= self.it_lim):
            print ("\tWarning:\n ST1 reached iteration limit.")
            
        if(np.pi*self.r_n_1**2 < self.fr.n_1*np.pi*self.r_in_orf_1**2):
            print("\tWarning:\n Outlet area smaller than total inlet area (dif = %.2f mm²)"%(self.fr.n_1*np.pi*self.r_in_orf_1**2 - np.pi*self.r_n_1**2))
            print("\tArea_in = %.2f mm²; \tArea_out = %.2f mm²"%(self.fr.n_1*np.pi*self.r_in_orf_1**2, np.pi*self.r_n_1**2))
        
        
        ##### ST1 design console report #####
        
        print("\n\nST1 Injector Geometry: \n")
        self.print_data("# of Inlet Channels", "%.1f"%(self.fr.n_1))
        self.print_data("Inlet Channel Radius", "%.3f \tmm"%(self.r_in_orf_1))
        self.print_data("Inlet Radius", "%.3f \tmm"%(self.r_in_pos_1))
        self.print_data("Swirl chamber Radius", "%.3f \tmm"%(self.r_s_1))
        self.print_data("Nozzle Radius", "%.3f \tmm"%(self.r_n_1))
        
        self.print_data("Nozzle Length", "%.3f \tmm"%(self.l_n_1))
        self.print_data("Swirl Chamber Length", "%.3f \tmm"%(self.l_s_1))
        self.print_data("Inlet Channel Length", "%.3f \tmm"%(self.l_in_1))
        
        print("\n\nST1 Injector Properties: ")
        self.print_data("Propellant", "%s"%(self.fr.propellant_1))
        self.print_data("Delta P", "%.2f bar"%(self.fr.delta_p_1))
        self.print_data("Mass flow target", "%.2f g/s"%(self.fr.m_1))
        
        self.print_data("A", "%.4f"%(self.A_1))
        self.print_data("A_eq", "%.4f"%(self.A_eq_1))
        self.print_data("K", "%.4f"%(self.K_1))
        self.print_data("Phi", "%.4f"%(self.phi_1))
        self.print_data("Phi_eq", "%.4f"%(self.phi_eq_1))
        self.print_data("a", "%.4f"%(self.a_1))
        self.print_data("a_eq", "%.4f"%(self.a_eq_1))
        self.print_data("Ksi", "%.4f"%(self.ksi_1))
        self.print_data("Lambda", "%.4f"%(self.lambda_1))
        self.print_data("mu_i", "%.4f"%(self.mu_i_1))
        self.print_data("C", "%.4f"%(self.C_1))
        #self.print_data("2Alpha ideal", "%.2f deg"%(self.alpha2_1))
        #self.print_data("2Alpha real", "%.2f deg"%(self.alpha2_eq_1))
        #self.print_data("2Alpha ideal jusante", "%.2f deg"%(self.alpha2_e_1))
        #self.print_data("2Alpha real jusante", "%.2f deg"%(self.alpha2_e_eq_1))
        self.print_data("2Alpha_avg","%.2f deg"%(self.alpha2_avg_1))
        
        self.print_data("r_mk_1", "%.2f mm"%(self.r_mk_1))
        self.print_data("r_mn_1", "%.2f mm"%(self.r_mn_1))
        self.print_data("t_fluid_1", "%.2f mm"%(self.t_fluid_1))
        self.print_data("Expected non-uniformity", "%.2f %%"%(self.I_1))


        ##### ST1 - Convergence plot #####
        
        self.step_x = np.ceil(self.contador/10)

        if(self.sh_plot_error == 1):
            fig, (ax1, ax2) = plt.subplots(2, dpi = self.dpi_value, figsize = [3.15, 3.15])
            fig.suptitle("Stage 1 - Solution convergence")
            
            ax1.plot(range(self.contador),self.inlet_diameter)
            ax1.set_xlabel("Iteration")
            ax1.set_ylabel("d_in [mm]")
            ax1.set_xticks(np.arange(0,self.contador+self.step_x,self.step_x))
            #ax1.grid()
            
            ax2.plot(range(self.contador),self.outlet_error)
            ax2.set_xlabel("Iteration")
            ax2.set_ylabel("R_n error [ad.]")
            ax2.set_xticks(np.arange(0,self.contador+self.step_x,self.step_x))
            #ax2.grid()

        
        ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
        ##### ##### ##### Parte 2 - Cálculo do ST2 #### ##### ##### ##### ##### #####
        ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
        
        if(self.fr.config == "Bi\n"):

            self.phi_2 = self.calc_phi(self.fr.alpha2_2)
            
            self.A_2 = self.calc_geom_char(self.phi_2)
            
            self.mu_2 = self.calc_mu(self.phi_2)
            
            #self.r_n_2 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(self.fr.m_2*1e-3/(self.mu_2*np.sqrt(2*self.fr.rho_2*self.fr.delta_p_2*1e5)))
            #self.r_n_2 = 1e3 * np.sqrt(self.fr.m_2*1e-3/(np.pi*self.mu_2*np.sqrt(2*self.fr.rho_2*self.fr.delta_p_2*1e5)))
            self.r_n_2 = 1e3*0.475*np.sqrt(self.fr.m_2*1e-3 / (self.mu_2 * np.sqrt(self.fr.rho_2*self.fr.delta_p_2*1e5)))
            
            print("\n\n\n\n\n ST2 Initial Values:")
            self.print_data("A_2", "%.2f"%(self.A_2))
            self.print_data("Phi_2", "%.3f"%(self.phi_2))
            self.print_data("R_n_2", "%.3f mm"%(self.r_n_2))
            self.print_data("2Alpha", "%.1f deg"%(self.fr.alpha2_2))
    
            self.inlet_diameter = []
            self.outlet_error = []
            
            
            self.erro_r = 100
            self.contador = 0
            
            if(self.show_it == 1): print("\n ST2 Iteration data:")
            while self.erro_r >= self.erro_max and self.contador < self.it_lim:
                
                
                ##### Iteração do injetor ideal #####
                
                self.r_in_pos_2 = self.fr.opening_coef_2*self.r_n_2
                
                self.r_in_orf_2 = np.sqrt((self.r_in_pos_2*self.r_n_2)/(self.fr.n_2*self.A_2))
                
                self.inlet_diameter.append(2*self.r_in_orf_2)
                
                self.r_s_2 = self.r_in_pos_2+self.r_in_orf_2
                
                self.l_in_2 = self.fr.ratio_l_in_2*2*self.r_in_orf_2
                
                self.l_n_2 = self.fr.ratio_l_n_2*2*self.r_n_2
    
                self.l_s_2 = self.fr.ratio_l_s_2*2*self.r_s_2
                
                self.a_2 = self.calc_a(self.phi_2)
                
                self.alpha2_2 = self.calc_2alpha(self.phi_2, "phi")
                
                self.alpha2_e_2 = self.calc_2alpha(self.a_2, "a")
                
                             
                
                ##### Iteração do injetor real #####
                
                self.rey_in_2 = self.calc_reynolds_in(self.fr.m_2, self.fr.din_visc_2, self.fr.n_2, self.r_in_orf_2, self.fr.rho_2)
                
                self.lambda_2 = self.calc_friction_coef(self.rey_in_2)
                
                self.ksi_2 = self.calc_hyd_loss(self.fr.in_type_2, self.r_s_2, self.l_in_2, self.fr.inlet_radius_2, self.r_in_orf_2, self.lambda_2)
                self.check_value("ksi_2", 0, float("inf"), self.ksi_2)
                
                self.A_eq_2 = self.calc_geom_char_eq(self.r_in_pos_2, self.r_n_2, self.fr.n_2, self.r_in_orf_2, self.lambda_2)
                
                self.phi_eq_2 = self.calc_phi_eq(self.A_eq_2)
                
                self.mu_eq_2 = self.calc_mu(self.phi_eq_2)
                                               
                self.mu_i_2 = self.calc_mu_i(self.mu_eq_2, self.A_2, self.r_in_pos_2, self.r_n_2, self.ksi_2)
                
                self.a_eq_2 = self.calc_a(self.phi_eq_2)
                
                self.alpha2_eq_2 = self.calc_2alpha(self.phi_eq_2, "phi")
                
                self.alpha2_e_eq_2 = self.calc_2alpha(self.a_eq_2, "a")
                
                
               
                ## cálculo do novo raio do bocal
                
                self.r_n_old = self.r_n_2
                
                #self.r_n_2 = 1e3/2* np.sqrt(4/np.pi) * np.sqrt(self.fr.m_2*1e-3/(self.mu_i_2*np.sqrt(2*self.fr.rho_2*self.fr.delta_p_2*1e5)))
                #self.r_n_2 = 1e3 * np.sqrt(self.fr.m_2*1e-3/(np.pi*self.mu_i_2*np.sqrt(2*self.fr.rho_2*self.fr.delta_p_2*1e5)))
                self.r_n_2 = 1e3*0.475*np.sqrt(self.fr.m_2*1e-3 / (self.mu_i_2 * np.sqrt(self.fr.rho_2*self.fr.delta_p_2*1e5)))
                
                self.erro_r = abs(self.r_n_old-self.r_n_2)/self.r_n_old
                
                self.outlet_error.append(self.erro_r)
                
                self.A_2 = self.r_in_pos_2*self.r_n_2/(self.fr.n_2*self.r_in_orf_2**2)
               
                
                self.phi_2 = self.calc_phi_eq(self.A_2)
                                                
                self.K_2 = self.A_eq_2/self.A_2
                self.check_value("K_2", 0, 1, self.K_2)
                
                self.r_mn_2 = self.r_n_2*np.sqrt(1-self.phi_eq_2)
                
                self.r_mk_2 = self.r_mn_2 * np.sqrt(2*(1-self.phi_eq_2)/(2-self.phi_eq_2))
                
                self.t_fluid_2 = self.r_n_2-self.r_mn_2
                
                self.C_2 = self.r_in_pos_2/self.r_n_2
                
                self.alpha2_avg_2 = self.calc_2alpha_avg(self.mu_i_2, self.A_eq_2, self.phi_eq_2, self.ksi_2, self.C_2)
                
                
                ## print dos dados da iteração
                if(self.show_it == 1):
                    print("\nIt. %d)"%(self.contador+1))
                    self.print_data("A_2", "%.4f"%(self.A_2))
                    self.print_data("A_eq_2", "%.4f"%(self.A_eq_2))
                    self.print_data("K_2", "%.4f"%(self.K_2))
                    self.print_data("Phi_2", "%.4f"%(self.phi_2))
                    self.print_data("Phi_eq_2", "%.4f"%(self.phi_eq_2))
                    self.print_data("a_2", "%.4f"%(self.a_2))
                    self.print_data("a_eq_2", "%.4f"%(self.a_eq_2))
                    self.print_data("rey_in_2", "%.4f"%(self.rey_in_2))
                    self.print_data("ksi_2", "%.4f"%(self.ksi_2))
                    self.print_data("lambda_2", "%.4f"%(self.lambda_2))
                    self.print_data("mu_i", "%.4f"%(self.mu_i_2))
                    self.print_data("C", "%.4f"%(self.C_2))
                    #self.print_data("2Alpha ideal", "%.2f deg"%(self.alpha2_2))
                    #self.print_data("2Alpha real", "%.2f deg"%(self.alpha2_eq_2))
                    #self.print_data("2Alpha ideal jusante", "%.2f deg"%(self.alpha2_e_2))
                    #self.print_data("2Alpha real jusante", "%.2f deg"%(self.alpha2_e_eq_2))
                    self.print_data("2Alpha_avg","%.2f deg"%(self.alpha2_avg_2))
                    self.print_data("R_n_2", "%.3f mm"%(self.r_n_2))
                    self.print_data("r_in_2", "%.3f mm"%(self.r_in_orf_2))
                    
                    self.print_data("Error", "%.2E"%(self.erro_r))
                
                self.contador += 1
            
            print("\n\n")
            
            ##Cálculo da não-uniformidade esperada do injetor
            self.I_2 = self.calc_non_uniformity(self.r_in_pos_2, self.r_n_2, self.fr.n_2, self.l_s_2, self.r_s_2, self.phi_eq_2)
                  
            
            if(self.contador == self.it_lim):
                print ("\tWarning:\n ST1 reached iteration limit.")
                
            if(np.pi*self.r_n_2**2 < self.fr.n_2*np.pi*self.r_in_orf_2**2):
                print("\tWarning:\n Outlet area smaller than total inlet area (dif = %.2f mm²)"%(self.fr.n_2*np.pi*self.r_in_orf_2**2 - np.pi*self.r_n_2**2))
                print("\tArea_in = %.2f mm²; \tArea_out = %.2f mm²"%(self.fr.n_2*np.pi*self.r_in_orf_2**2, np.pi*self.r_n_2**2))
            
            
            ##### ST2 design console report #####
            
            print("\n\nST2 Injector Geometry: \n")
            self.print_data("# of Inlet Channels", "%.1f"%(self.fr.n_2))
            self.print_data("Inlet Channel Radius", "%.3f \tmm"%(self.r_in_orf_2))
            self.print_data("Inlet Radius", "%.3f \tmm"%(self.r_in_pos_2))
            self.print_data("Swirl chamber Radius", "%.3f \tmm"%(self.r_s_2))
            self.print_data("Nozzle Radius", "%.3f \tmm"%(self.r_n_2))
            
            self.print_data("Nozzle Length", "%.3f \tmm"%(self.l_n_2))
            self.print_data("Swirl Chamber Length", "%.3f \tmm"%(self.l_s_2))
            self.print_data("Inlet Channel Length", "%.3f \tmm"%(self.l_in_2))
            
            print("\n\nST2 Injector Properties: ")
            self.print_data("Propellant", "%s"%(self.fr.propellant_2))
            self.print_data("Delta P", "%.2f bar"%(self.fr.delta_p_2))
            self.print_data("Mass flow target", "%.2f g/s"%(self.fr.m_2))
            
            self.print_data("A", "%.4f"%(self.A_2))
            self.print_data("A_eq", "%.4f"%(self.A_eq_2))
            self.print_data("K", "%.4f"%(self.K_2))
            self.print_data("Phi", "%.4f"%(self.phi_2))
            self.print_data("Phi_eq", "%.4f"%(self.phi_eq_2))
            self.print_data("a", "%.4f"%(self.a_2))
            self.print_data("a_eq", "%.4f"%(self.a_eq_2))
            self.print_data("Ksi", "%.4f"%(self.ksi_2))
            self.print_data("Lambda", "%.4f"%(self.lambda_2))
            self.print_data("mu_i", "%.4f"%(self.mu_i_2))
            self.print_data("C", "%.4f"%(self.C_2))
            #self.print_data("2Alpha ideal", "%.2f deg"%(self.alpha2_2))
            #self.print_data("2Alpha real", "%.2f deg"%(self.alpha2_eq_2))
            #self.print_data("2Alpha ideal jusante", "%.2f deg"%(self.alpha2_e_2))
            #self.print_data("2Alpha real jusante", "%.2f deg"%(self.alpha2_e_eq_2))
            self.print_data("2Alpha_avg","%.2f deg"%(self.alpha2_avg_2))
            
            self.print_data("r_mk_2", "%.2f mm"%(self.r_mk_2))
            self.print_data("r_mn_2", "%.2f mm"%(self.r_mn_2))
            self.print_data("t_fluid_2", "%.2f mm"%(self.t_fluid_2))
            self.print_data("Expected non-uniformity", "%.2f %%"%(self.I_2))
            
            
            ##### ST2 - Convergence plot #####
            
            self.step_x = np.ceil(self.contador/10)
            
            if(self.sh_plot_error == 1):
                fig, (ax3, ax4) = plt.subplots(2, dpi = self.dpi_value)
                fig.suptitle("Stage 2 - Solution convergence")
                
                ax3.plot(range(self.contador),self.inlet_diameter)
                ax3.set_xlabel("Iteration")
                ax3.set_ylabel("d_in [mm]")
                ax3.set_xticks(np.arange(0,self.contador+self.step_x,self.step_x))
                #ax3.grid()
                
                ax4.plot(range(self.contador),self.outlet_error)
                ax4.set_xlabel("Iteration")
                ax4.set_ylabel("R_n error [ad.]")
                ax4.set_xticks(np.arange(0,self.contador+self.step_x,self.step_x))
                #ax4.grid()
            
        
            ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
            ##### ##### ##### Parte 3 - Integração dos Estágios # ##### ##### ##### #####
            ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
            
            self.gap = self.r_mk_2 - (self.r_n_1 + self.fr.t_w)
            
            self.fluid_t_2_s = self.r_s_2 - self.r_mk_2
            
            print("\n\nGeometry Check: \t")
            self.print_data("D_n_1", "%.3f mm"%(2*self.r_n_1))
            self.print_data("D_w_1", "%.3f mm"%(2*self.r_n_1+2*self.fr.t_w))
            self.print_data("D_n_2", "%.3f mm"%(2*self.r_n_2))
            self.print_data("fluid_t_2_s", "%.3f mm"%(self.fluid_t_2_s))
            
            
            if(self.gap >= self.fr.delta/100*self.fluid_t_2_s):
                print("\n\tHydraulic Independance OK! \n\tDeltaR = %.3f mm (mín: %.2f mm)"%(self.gap, self.fr.delta/100*self.fluid_t_2_s))
                
            else:
                print("\n\tERROR: Hydraulic Independance FAIL! \n\tDeltaR = %.3f mm (mín: %.2f mm)"%(self.gap, self.fr.delta/100*self.fluid_t_2_s))
                self.check_out = 0
                           
                    ##### ##### ##### ##### ######
                    # Alteração do ângulo do ST2 #
                    ##### ##### ##### ##### ######
                    
            if(self.fr.angle_dif == "disabled"):
                self.check_out = -1
                
            else:
                    
                    #self.alpha2_goal = 2*(self.alpha2_eq_1/2 - self.fr.angle_dif)
                    self.alpha2_goal = 2*(self.alpha2_e_eq_1/2 - self.fr.angle_dif)
                    
                    #if(self.alpha2_goal < self.alpha2_eq_2):
                    if(self.alpha2_goal < self.alpha2_e_eq_2):
                    
                        #self.a_2 = 2*(1-self.phi_eq_2)**2/(2-self.phi_eq_2)
                        #self.a_2 = 2*(1-self.phi_2)**2/(2-self.phi_2)
                        
                        self.r_out = self.r_n_2*np.sqrt(self.a_2*(1+1/(np.tan(self.alpha2_goal/2*np.pi/180))**2))
                        
                        self.l_out = self.fr.ratio_l_n_out*2*self.r_n_2
                        
                        #self.phi_out = self.phi_eq_2/np.sqrt(3-2*self.phi_eq_2)
                        self.phi_out = self.phi_2/np.sqrt(3-2*self.phi_2)
                        
                        self.r_mn_out = self.r_out*np.sqrt(1-self.phi_out)
                        
                        self.t_fluid_out = self.r_out-self.r_mn_out
                        
                        
                        self.u_tot_2 = (self.fr.m_2*1e-3)/(self.fr.rho_2*self.fr.n_2*np.pi*(self.r_in_orf_2*1e-3)**2)
                        
                        self.u_wave_2 = np.sqrt(1-self.a_2) * self.u_tot_2
                        
                        self.u_a_out = np.sqrt(1-self.a_2/(self.r_out/self.r_n_2)**2) * self.u_tot_2
                        
                        #u_th = u_an = "axial velocity in the nozzle"
                        self.u_an_2 = np.sqrt(self.phi_eq_2/(2-self.phi_eq_2)) * self.u_tot_2
                        
                        
                        self.lambda_i = self.u_a_out/self.u_an_2

                        self.lambda_i_max = np.sqrt((2-self.phi_2)/self.phi_2)
                        
                        if(self.lambda_i > self.lambda_i_max):
                            print("Error, Lambda exceeds maximum allowable value (Outlet radius too large)!")
                            self.print_data("lambda_i", "%.3f"%(self.lambda_i))
                            self.print_data("lambda_i_max", "%.3f"%(self.lambda_i_max))
                            
                        
                        
                        if(self.r_out <= self.r_n_2):
                            
                            print("\n\tError: Invalid ST2 exit radius.")
                            
                            self.check_out = -1
                            
                        else:
                            print("\n\nST2 Modified Injector Geometry: \t")
                            self.print_data("2Alpha_2_goal", "%.1f deg"%(self.alpha2_goal))
                            self.print_data("Phi_out", "%.3f"%(self.phi_out))
                            self.print_data("R_out", "%.3f mm"%(self.r_out))
                            self.print_data("L_out", "%.3f mm"%(self.l_out))
                            self.print_data("R_mn_out", "%.3f mm"%(self.r_mn_out))
                            self.print_data("t_fluid_out", "%.3f mm"%(self.t_fluid_out))
                            self.print_data("u_tot_2", "%.3f m/s"%(self.u_tot_2))
                            self.print_data("u_wave_2", "%.3f m/s"%(self.u_wave_2))
                            self.print_data("u_a_out", "%.3f m/s"%(self.u_a_out))
                            self.print_data("u_an_2", "%.3f m/s"%(self.u_an_2))
                            self.print_data("lambda_i", "%.3f"%(self.lambda_i))
                            self.print_data("lambda_i_max", "%.3f"%(self.lambda_i_max))
                            
                            
                            self.check_out = 1
                        
                    else:
                        print("\n\tWarning: Intercepting angle greater than input. \n\tInt_angle = %.1f deg"%(self.alpha2_eq_1/2 - self.alpha2_eq_2/2))
                        
                        self.check_out = -1

        


        ##### ##### ##### ##### #
        ##### Injector plot #####
        ##### ##### ##### ##### #
        
        if(self.fr.config == "Mono\n"):
        
            self.alpha2_plot_1 = self.alpha2_avg_1

            ##### conventional injector plot #####
            
            ## ST1

            self.st1_n_bottom = 0
            self.st1_n_top = self.st1_n_bottom + self.l_n_1
            self.st1_s_bottom = self.st1_n_top + (self.r_s_1-self.r_n_1)/np.tan(self.fr.trans_angle_1*np.pi/180)
            self.st1_s_top = self.st1_s_bottom + self.l_s_1
 
            
            self.st1_x = [self.r_n_1, self.r_n_1, self.r_s_1, self.r_s_1, 0]          
            self.st1_y = [self.st1_n_bottom, self.st1_n_top, self.st1_s_bottom, self.st1_s_top, self.st1_s_top]
            
            ## pontos do orifício tangencial do ST1
            
            self.st1_in_top = []
            self.st1_in_bottom = []
            
            self.circ_x_1 = np.arange(self.r_s_1-2*self.r_in_orf_1,self.r_s_1,0.005)
            for x in self.circ_x_1:
                self.st1_in_top.append((self.st1_s_top-self.r_in_orf_1) + np.sqrt(abs(self.r_in_orf_1**2-(x-(self.r_s_1-self.r_in_orf_1))**2)))
                self.st1_in_bottom.append((self.st1_s_top-self.r_in_orf_1) - np.sqrt(abs(self.r_in_orf_1**2-(x-(self.r_s_1-self.r_in_orf_1))**2)))
            
            ## pontos do fluido interno do st1
            self.st1_fluid_n_bottom = 0
            self.st1_fluid_n_top = self.st1_n_bottom + self.l_n_1
            self.st1_fluid_s_bottom = self.st1_n_top + (self.r_s_1-self.r_n_1)/np.tan(self.fr.trans_angle_1*np.pi/180)
            self.st1_fluid_s_top = self.st1_s_bottom + self.l_s_1            
            
            self.st1_fluid_x = [self.r_mn_1, self.r_mn_1, self.r_mk_1, self.r_mk_1]
            self.st1_fluid_y = [self.st1_fluid_n_bottom, self.st1_fluid_n_top, self.st1_fluid_s_bottom, self.st1_fluid_s_top]
                        
            
            ## pontos do spray do st1
            
            self.st1_spray_y_i = 0
            self.st1_spray_y_f = -1*self.st1_s_top*1.1
            self.st1_spray_y_step = (self.st1_spray_y_f-self.st1_spray_y_i)/100
            self.st1_spray_y = np.arange(self.st1_spray_y_i,self.st1_spray_y_f+self.st1_spray_y_step,self.st1_spray_y_step)
            
            self.st1_spray_x = []
            
            self.r_av_1 = (self.r_n_1 + self.r_mn_1)/2         
                    
            for y in self.st1_spray_y:
                self.st1_spray_x.append(self.r_av_1*np.sqrt(1+y**2/(self.r_av_1*np.tan((90-0.5*self.alpha2_plot_1)*np.pi/180))**2))
            
            if(self.sh_plot_inj == 1):
                
                fig, ax5 = plt.subplots(dpi = self.dpi_value)
                ax5.plot(self.st1_x,self.st1_y, linewidth = 2, color = '0.3', linestyle = '-')
                ax5.plot(self.circ_x_1,self.st1_in_top, linewidth = 1, color = '0.3', linestyle = '-')
                ax5.plot(self.circ_x_1,self.st1_in_bottom, linewidth = 1, color = '0.3', linestyle = '-')
                ax5.plot(self.st1_fluid_x, self.st1_fluid_y, linewidth = 1, color = 'b', linestyle = '--')
                ax5.plot(self.st1_spray_x, self.st1_spray_y, linewidth = 2, color = 'b', linestyle = '--')
                #ax5.grid()
                ax5.set_xlabel("[mm]")
                ax5.set_ylabel("[mm]")
                
                self.x_min = 0
                self.x_max = 2*self.st1_s_top*1.1
                self.y_min = -1*self.st1_s_top*1.1
                self.y_max = self.st1_s_top*1.1
                
                ax5.set_xlim(self.x_min, self.x_max)
                ax5.set_ylim(self.y_min, self.y_max)
                ax5.set_aspect("equal")
                
                          
                self.st1_legend = "ST1:\n   R_n: %.2f mm\n   R_in: %.2f mm\n   R_s: %.2f mm\n   r_in: %.2f mm\n   n_in: %d \n"%(self.r_n_1, self.r_in_pos_1, self.r_s_1, self.r_in_orf_1, self.fr.n_1)
                self.bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
                ax5.text(self.r_s_1*1.5, 0, self.st1_legend, fontsize=self.fontsize, ha="left", va="bottom", bbox=self.bbox_props)

                
        
        ##### plots para injetores bipropelente #####
        
        
        elif(self.fr.config == "Bi\n"):
            
            self.alpha2_plot_1 = self.alpha2_avg_1
            self.alpha2_plot_2 = self.alpha2_avg_2
            
            if(self.sh_plot_inj == 1 and (self.check_out == 1 or self.check_out == -1)):
                ##### conventional injector plot #####
                
                ## ST1
                
                self.st1_n_bottom = self.fr.recess
                self.st1_n_top = self.st1_n_bottom + self.l_n_1
                self.st1_s_bottom = self.st1_n_top + (self.r_s_1-self.r_n_1)/np.tan(self.fr.trans_angle_1*np.pi/180)
                self.st1_s_top = self.st1_s_bottom + self.l_s_1
                
                self.st1_x = [self.r_n_1, self.r_n_1, self.r_s_1, self.r_s_1, 0]
                self.st1_y = [self.st1_n_bottom, self.st1_n_top, self.st1_s_bottom, self.st1_s_top, self.st1_s_top]
                
                
                ## pontos do orifício tangencial do ST1
                
                self.st1_in_top = []
                self.st1_in_bottom = []
                
                self.circ_x_1 = np.arange(self.r_s_1-2*self.r_in_orf_1,self.r_s_1,0.005)
                for x in self.circ_x_1:
                    self.st1_in_top.append((self.st1_s_top-self.r_in_orf_1) + np.sqrt(abs(self.r_in_orf_1**2-(x-(self.r_s_1-self.r_in_orf_1))**2)))
                    self.st1_in_bottom.append((self.st1_s_top-self.r_in_orf_1) - np.sqrt(abs(self.r_in_orf_1**2-(x-(self.r_s_1-self.r_in_orf_1))**2)))
    
                
                ## pontos do fluido interno do st1
                self.st1_fluid_n_bottom = self.fr.recess
                self.st1_fluid_n_top = self.st1_n_bottom + self.l_n_1
                self.st1_fluid_s_bottom = self.st1_n_top + (self.r_s_1-self.r_n_1)/np.tan(self.fr.trans_angle_1*np.pi/180)
                self.st1_fluid_s_top = self.st1_s_bottom + self.l_s_1            
                
                self.st1_fluid_x = [self.r_mn_1, self.r_mn_1, self.r_mk_1, self.r_mk_1]
                self.st1_fluid_y = [self.st1_fluid_n_bottom, self.st1_fluid_n_top, self.st1_fluid_s_bottom, self.st1_fluid_s_top]
                
                ## pontos do spray do st1            
                self.st1_spray_y_i = 0
                self.st1_spray_y_f = -1*self.st1_s_top*1.1 - self.fr.recess
                self.st1_spray_y_step = (self.st1_spray_y_f-self.st1_spray_y_i)/100
                self.st1_spray_y = np.arange(self.st1_spray_y_i,self.st1_spray_y_f+self.st1_spray_y_step,self.st1_spray_y_step)
                
                self.st1_spray_x = []
                
                self.r_av_1 = (self.r_n_1 + self.r_mn_1)/2 
                
                
                        
                for y in self.st1_spray_y:
                    self.st1_spray_x.append(self.r_av_1*np.sqrt(1+y**2/(self.r_av_1*np.tan((90-0.5*self.alpha2_plot_1)*np.pi/180))**2))
                    
                for i in range(len(self.st1_spray_y)):
                    self.st1_spray_y[i] = self.st1_spray_y[i] + self.fr.recess
                
                
                
                
                ## ST2
                            
                self.st2_n_bottom = 0
                self.st2_n_top = self.st2_n_bottom + self.l_n_2
                self.st2_s_bottom = self.st2_n_top + (self.r_s_2-self.r_n_2)/np.tan(self.fr.trans_angle_2*np.pi/180)
                self.st2_s_top = self.st2_s_bottom + self.l_s_2
                
                            
                ## pontos da estrutura do ST2
                self.st2_x = [self.r_n_2, self.r_n_2, self.r_s_2, self.r_s_2, self.r_n_1+self.fr.t_w, self.r_n_1+self.fr.t_w, self.r_n_1]
                self.st2_y = [self.st2_n_bottom, self.st2_n_top, self.st2_s_bottom, self.st2_s_top, self.st2_s_top, self.fr.recess, self.fr.recess]
                
                
                ## pontos do orifício tangencial do ST2
                
                self.st2_in_top = []
                self.st2_in_bottom = []
                
                self.circ_x_2 = np.arange(self.r_s_2-2*self.r_in_orf_2,self.r_s_2,0.005)
                for x in self.circ_x_2:
                    self.st2_in_top.append((self.st2_s_top-self.r_in_orf_2) + np.sqrt(abs(self.r_in_orf_2**2-(x-(self.r_s_2-self.r_in_orf_2))**2)))
                    self.st2_in_bottom.append((self.st2_s_top-self.r_in_orf_2) - np.sqrt(abs(self.r_in_orf_2**2-(x-(self.r_s_2-self.r_in_orf_2))**2)))
                    
                ## pontos do fluido interno do ST2
                self.st2_fluid_n_bottom = 0
                self.st2_fluid_n_top = self.st2_fluid_n_bottom + self.l_n_2
                self.st2_fluid_s_bottom = self.st2_fluid_n_top + (self.r_s_2-self.r_n_2)/np.tan(self.fr.trans_angle_2*np.pi/180)
                self.st2_fluid_s_top = self.st2_fluid_s_bottom + self.l_s_2            
                
                self.st2_fluid_x = [self.r_mn_2, self.r_mn_2, self.r_mk_2, self.r_mk_2]
                self.st2_fluid_y = [self.st2_fluid_n_bottom, self.st2_fluid_n_top, self.st2_fluid_s_bottom, self.st2_fluid_s_top]
                            
                
                ## pontos do spray do ST2
                
                self.st2_spray_y_i = 0
                self.st2_spray_y_f = -1*self.st1_s_top*1.1
                self.st2_spray_y_step = (self.st2_spray_y_f-self.st2_spray_y_i)/100
                self.st2_spray_y = np.arange(self.st2_spray_y_i,self.st2_spray_y_f+self.st2_spray_y_step,self.st2_spray_y_step)
                
                self.st2_spray_x = []
                
                self.r_av_2 = (self.r_n_2 + self.r_mn_2)/2         
                        
                for y in self.st2_spray_y:
                    self.st2_spray_x.append(self.r_av_2*np.sqrt(1+y**2/(self.r_av_2*np.tan((90-0.5*self.alpha2_plot_2)*np.pi/180))**2))
                    
                
                
                
    
    
                fig, ax6 = plt.subplots(dpi = self.dpi_value)
                
                ax6.plot(self.st1_x,self.st1_y, linewidth = 1, color = '0.3', linestyle = '-')
                ax6.plot(self.circ_x_1,self.st1_in_top, linewidth = 1, color = '0.3', linestyle = '-')
                ax6.plot(self.circ_x_1,self.st1_in_bottom, linewidth = 1, color = '0.3', linestyle = '-')
                ax6.plot(self.st1_fluid_x, self.st1_fluid_y, linewidth = 1, color = 'b', linestyle = '--')
                ax6.plot(self.st1_spray_x, self.st1_spray_y, linewidth = 2, color = 'b', linestyle = '--')
                
                ax6.plot(self.st2_x,self.st2_y, linewidth = 1, color = '0.3', linestyle = '-')
                ax6.plot(self.circ_x_2,self.st2_in_top, linewidth = 1, color = '0.3', linestyle = '-')
                ax6.plot(self.circ_x_2,self.st2_in_bottom, linewidth = 1, color = '0.3', linestyle = '-')
                ax6.plot(self.st2_fluid_x, self.st2_fluid_y, linewidth = 1, color = 'y', linestyle = '--')
                ax6.plot(self.st2_spray_x, self.st2_spray_y, linewidth = 2, color = 'y', linestyle = '--')
                
                #ax6.grid()
                ax6.set_xlabel("[mm]")
                ax6.set_ylabel("[mm]")
                
                self.x_min = 0
                self.x_max = 2*self.st1_s_top*1.1
                self.y_min = -1*self.st1_s_top*1.1
                self.y_max = self.st1_s_top*1.1
                
                ax6.set_xlim(self.x_min, self.x_max)
                ax6.set_ylim(self.y_min, self.y_max)
                ax6.set_aspect("equal")
                
                              
                #self.st1_legend = "ST1:\n   R_n: %.2f mm\n   R_in: %.2f mm\n   R_s: %.2f mm\n   r_in: %.2f mm\n   n_in: %d \n   recess: %.2f mm"%(self.r_n_1, self.r_in_pos_1, self.r_s_1, self.r_in_orf_1, self.fr.n_1, self.fr.recess)
                #self.st2_legend = "ST2:\n   R_n: %.2f mm\n   R_in: %.2f mm\n   R_s: %.2f mm\n   r_in: %.2f mm\n   n_in: %d \n   int_angle: %.1f°"%(self.r_n_2, self.r_in_pos_2, self.r_s_2, self.r_in_orf_2, self.fr.n_2, self.alpha2_eq_1/2 - self.alpha2_eq_2/2)
                #self.bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
                #ax6.text(max([self.r_s_1,self.r_s_2])*1.5, 0, self.st1_legend, fontsize=self.fontsize, ha="left", va="bottom", bbox=self.bbox_props)
                #ax6.text(max([self.r_s_1,self.r_s_2])*1.5+self.offset, 0, self.st2_legend, fontsize=self.fontsize, ha="left", va="bottom", bbox=self.bbox_props)
            
        
        
            
                if(self.check_out == 1):
                    ##### modified injector plot #####
                    
                    ## ST1
                    
                                
                    self.st1_n_bottom = self.fr.recess
                    self.st1_n_top = self.st1_n_bottom + self.l_n_1
                    self.st1_s_bottom = self.st1_n_top + (self.r_s_1-self.r_n_1)/np.tan(self.fr.trans_angle_1*np.pi/180)
                    self.st1_s_top = self.st1_s_bottom + self.l_s_1
                    
                    self.st1_x = [self.r_n_1, self.r_n_1, self.r_s_1, self.r_s_1, 0]
                    self.st1_y = [self.st1_n_bottom, self.st1_n_top, self.st1_s_bottom, self.st1_s_top, self.st1_s_top]
                    
                    
                    ## pontos do orifício tangencial do ST1
                    
                    self.st1_in_top = []
                    self.st1_in_bottom = []
                    
                    self.circ_x_1 = np.arange(self.r_s_1-2*self.r_in_orf_1,self.r_s_1,0.005)
                    for x in self.circ_x_1:
                        self.st1_in_top.append((self.st1_s_top-self.r_in_orf_1) + np.sqrt(abs(self.r_in_orf_1**2-(x-(self.r_s_1-self.r_in_orf_1))**2)))
                        self.st1_in_bottom.append((self.st1_s_top-self.r_in_orf_1) - np.sqrt(abs(self.r_in_orf_1**2-(x-(self.r_s_1-self.r_in_orf_1))**2)))
                    
                    ## pontos do fluido interno do st1
                    self.st1_fluid_n_bottom = self.fr.recess
                    self.st1_fluid_n_top = self.st1_n_bottom + self.l_n_1
                    self.st1_fluid_s_bottom = self.st1_n_top + (self.r_s_1-self.r_n_1)/np.tan(self.fr.trans_angle_1*np.pi/180)
                    self.st1_fluid_s_top = self.st1_s_bottom + self.l_s_1            
                    
                    self.st1_fluid_x = [self.r_mn_1, self.r_mn_1, self.r_mk_1, self.r_mk_1]
                    self.st1_fluid_y = [self.st1_fluid_n_bottom, self.st1_fluid_n_top, self.st1_fluid_s_bottom, self.st1_fluid_s_top]
                    
                    ## pontos do spray do st1            
                    self.st1_spray_y_i = 0
                    self.st1_spray_y_f = -1*self.st1_s_top*1.1 - self.fr.recess
                    self.st1_spray_y_step = (self.st1_spray_y_f-self.st1_spray_y_i)/100
                    self.st1_spray_y = np.arange(self.st1_spray_y_i,self.st1_spray_y_f+self.st1_spray_y_step,self.st1_spray_y_step)
                    
                    self.st1_spray_x = []
                    
                    self.r_av_1 = (self.r_n_1 + self.r_mn_1)/2         
                            
                    for y in self.st1_spray_y:
                        self.st1_spray_x.append(self.r_av_1*np.sqrt(1+y**2/(self.r_av_1*np.tan((90-0.5*self.alpha2_plot_1)*np.pi/180))**2))
                        
                    for i in range(len(self.st1_spray_y)):
                        self.st1_spray_y[i] = self.st1_spray_y[i] + self.fr.recess
                    
                    
                    
                    
                    ## ST2
                    
                    self.st2_out_bottom = 0
                    self.st2_out_top = self.l_out
                    self.st2_n_bottom = self.st2_out_top + (self.r_out-self.r_n_2)/np.tan(self.fr.trans_angle_out*np.pi/180)
                    self.st2_n_top = self.st2_n_bottom + self.l_n_2
                    self.st2_s_bottom = self.st2_n_top + (self.r_s_2-self.r_n_2)/np.tan(self.fr.trans_angle_2*np.pi/180)
                    self.st2_s_top = self.st2_s_bottom + self.l_s_2
                                
                    ## pontos da estrutura do ST2
                    self.st2_x = [self.r_out, self.r_out, self.r_n_2, self.r_n_2, self.r_s_2, self.r_s_2, self.r_n_1+self.fr.t_w, self.r_n_1+self.fr.t_w, self.r_n_1]
                    self.st2_y = [self.st2_out_bottom, self.st2_out_top, self.st2_n_bottom, self.st2_n_top, self.st2_s_bottom, self.st2_s_top, self.st2_s_top, self.fr.recess, self.fr.recess]
                    
                    
                    ## pontos do orifício tangencial do ST2
                    
                    self.st2_in_top = []
                    self.st2_in_bottom = []
                    
                    self.circ_x_2 = np.arange(self.r_s_2-2*self.r_in_orf_2,self.r_s_2,0.005)
                    for x in self.circ_x_2:
                        self.st2_in_top.append((self.st2_s_top-self.r_in_orf_2) + np.sqrt(abs(self.r_in_orf_2**2-(x-(self.r_s_2-self.r_in_orf_2))**2)))
                        self.st2_in_bottom.append((self.st2_s_top-self.r_in_orf_2) - np.sqrt(abs(self.r_in_orf_2**2-(x-(self.r_s_2-self.r_in_orf_2))**2)))
        
                    
                    ## pontos do fluido interno do ST2
                    self.st2_fluid_out_bottom = 0
                    self.st2_fluid_out_top = self.l_out
                    self.st2_fluid_n_bottom = self.st2_out_top + (self.r_out-self.r_n_2)/np.tan(self.fr.trans_angle_out*np.pi/180)
                    self.st2_fluid_n_top = self.st2_n_bottom + self.l_n_2
                    self.st2_fluid_s_bottom = self.st2_n_top + (self.r_s_2-self.r_n_2)/np.tan(self.fr.trans_angle_2*np.pi/180)
                    self.st2_fluid_s_top = self.st2_s_bottom + self.l_s_2            
                    
                    self.st2_fluid_x = [self.r_mn_out, self.r_mn_out, self.r_mn_2, self.r_mn_2, self.r_mk_2, self.r_mk_2]
                    self.st2_fluid_y = [self.st2_fluid_out_bottom, self.st2_fluid_out_top, self.st2_fluid_n_bottom, self.st2_fluid_n_top, self.st2_fluid_s_bottom, self.st2_fluid_s_top]
                                
                    
                    ## pontos do spray do ST2
                    
                    self.st2_spray_y_i = 0
                    self.st2_spray_y_f = -1*self.st1_s_top*1.1
                    self.st2_spray_y_step = (self.st2_spray_y_f-self.st2_spray_y_i)/100
                    self.st2_spray_y = np.arange(self.st2_spray_y_i,self.st2_spray_y_f+self.st2_spray_y_step,self.st2_spray_y_step)
                    
                    self.st2_spray_x = []
                    
                    self.r_av_2 = (self.r_out + self.r_mn_out)/2    
                            
                    for y in self.st2_spray_y:
                        self.st2_spray_x.append(self.r_av_2*np.sqrt(1+y**2/(self.r_av_2*np.tan((90-0.5*self.alpha2_goal)*np.pi/180))**2))
                    
                    
                    fig, ax7 = plt.subplots(dpi = self.dpi_value)
                    
                    ax7.plot(self.st1_x,self.st1_y, linewidth = 1, color = '0.3', linestyle = '-')
                    ax7.plot(self.circ_x_1,self.st1_in_top, linewidth = 1, color = '0.3', linestyle = '-')
                    ax7.plot(self.circ_x_1,self.st1_in_bottom, linewidth = 1, color = '0.3', linestyle = '-')
                    ax7.plot(self.st1_fluid_x, self.st1_fluid_y, linewidth = 1, color = 'b', linestyle = '--')
                    ax7.plot(self.st1_spray_x, self.st1_spray_y, linewidth = 2, color = 'b', linestyle = '--')
                    
                    ax7.plot(self.st2_x,self.st2_y, linewidth = 1, color = '0.3', linestyle = '-')
                    ax7.plot(self.circ_x_2,self.st2_in_top, linewidth = 1, color = '0.3', linestyle = '-')
                    ax7.plot(self.circ_x_2,self.st2_in_bottom, linewidth = 1, color = '0.3', linestyle = '-')
                    ax7.plot(self.st2_fluid_x, self.st2_fluid_y, linewidth = 1, color = 'y', linestyle = '--')
                    ax7.plot(self.st2_spray_x, self.st2_spray_y, linewidth = 2, color = 'y', linestyle = '--')
                    
                    #ax7.grid()
                    ax7.set_xlabel("[mm]")
                    ax7.set_ylabel("[mm]")
                    
                    self.x_min = 0
                    self.x_max = 2*self.st1_s_top*1.1
                    self.y_min = -1*self.st1_s_top*1.1
                    self.y_max = self.st1_s_top*1.1
                    
                    ax7.set_xlim(self.x_min, self.x_max)
                    ax7.set_ylim(self.y_min, self.y_max)
                    ax7.set_aspect("equal")
                      
                    #self.st1_legend = "ST1:\n   R_n: %.2f mm\n   R_in: %.2f mm\n   R_s: %.2f mm\n   r_in: %.2f mm\n   n_in: %d \n   recess: %.2f mm"%(self.r_n_1, self.r_in_pos_1, self.r_s_1, self.r_in_orf_1, self.fr.n_1, self.fr.recess)
                    #self.st2_legend = "ST2:\n   R_n: %.2f mm\n   R_in: %.2f mm\n   R_s: %.2f mm\n   r_in: %.2f mm\n   n_in: %d \n   R_out: %.2f mm"%(self.r_n_2, self.r_in_pos_2, self.r_s_2, self.r_in_orf_2, self.fr.n_2, self.r_out)
                    #self.bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
                    #ax7.text(max([self.r_s_1,self.r_s_2])*1.5, 0, self.st1_legend, fontsize=self.fontsize, ha="left", va="bottom", bbox=self.bbox_props)
                    #ax7.text(max([self.r_s_1,self.r_s_2])*1.5+self.offset, 0, self.st2_legend, fontsize=self.fontsize, ha="left", va="bottom", bbox=self.bbox_props)
                
                
                elif(self.check_out == -1):
                    print("\n\nPlot: Modified injector not plotted\n\t(Angle modification not needed or disabled by user)")

            elif(self.check_out == 0):
                print("\n\nPlot: injector not plotted\n\t(Hydraulic Independance FAIL!)")
                
            elif(self.sh_plot_inj == 0):
                print("\n\nPlot: injector not plotted\n\t(Disabled by user)")
                
            else:
                print("\n\nPlot: injector not plotted\n\t(Unknow Error)")
            
        elif(self.sh_plot_inj == 1):
            print("Plot: configuration not recognized")

        plt.show()