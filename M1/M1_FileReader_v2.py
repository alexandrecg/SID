# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 2020

@author: A. Goulart
"""

class FileReader_1:
    def setFolderName(self,fo_name):
        self.folder_name = fo_name
    def setFileName(self,fi_name):
        self.file_name = fi_name
    def read(self):

        ##### Injector parameters reading #####
        
        file_loc = str(self.folder_name) + '/' + str(self.file_name) + '.txt'
        
        Input = open(str(file_loc),'r') 
        
        for i in range(5): Input.readline()
        
        injector_ID = str(Input.readline())
        print("Identificação do injetor: \t", injector_ID)
        
        Input.readline()
        self.config = str(Input.readline())
        print("Configuração do injetor: \t", self.config)
        
        
        
        if(self.config == "Mono\n"):
           
        
        ##### Beginning of ST1 data if monopropellant #####
        
            Input.readline()    
            self.in_type_1 = str(Input.readline())
            
            Input.readline()
            self.inlet_radius_1 = float(Input.readline())
            
            Input.readline()    
            self.propellant_1 = str(Input.readline())
            
            Input.readline()    
            self.m_1 = float(Input.readline())
            
            Input.readline()    
            self.delta_p_1 = float(Input.readline())
            
            Input.readline()    
            self.alpha2_1 = float(Input.readline())
            
            Input.readline()    
            self.n_1 = int(Input.readline())
            
            Input.readline()    
            self.opening_coef_1 = float(Input.readline())
            
            Input.readline()
            self.trans_angle_1 = float(Input.readline())
            
            Input.readline()    
            self.ratio_l_in_1 = float(Input.readline())
            
            Input.readline()    
            self.ratio_l_s_1 = float(Input.readline())
            
            Input.readline()    
            self.ratio_l_n_1 = float(Input.readline())
        
            print("ST1: Data read Successful")        
        
        elif(self.config == "Bi\n"):
            
            Input.readline()
            self.t_w = float(Input.readline())
            
            Input.readline()
            self.delta = float(Input.readline())
            
            Input.readline()
            self.angle_dif = float(Input.readline())
            
            Input.readline()
            self.recess = float(Input.readline())
            
            ##### Beginning of ST1 data if bipropellant #####
        
            Input.readline()    
            self.in_type_1 = str(Input.readline())
            
            Input.readline()
            self.inlet_radius_1 = float(Input.readline())
            
            Input.readline()    
            self.propellant_1 = str(Input.readline())
            
            Input.readline()    
            self.m_1 = float(Input.readline())
            
            Input.readline()    
            self.delta_p_1 = float(Input.readline())
            
            Input.readline()    
            self.alpha2_1 = float(Input.readline())
            
            Input.readline()    
            self.n_1 = int(Input.readline())
            
            Input.readline()    
            self.opening_coef_1 = float(Input.readline())
            
            Input.readline()
            self.trans_angle_1 = float(Input.readline())
            
            Input.readline()    
            self.ratio_l_in_1 = float(Input.readline())
            
            Input.readline()    
            self.ratio_l_s_1 = float(Input.readline())
            
            Input.readline()    
            self.ratio_l_n_1 = float(Input.readline())
            
            print("ST1: Data read Successful")   
            
            ##### Beginning of ST2 data #####
            
            Input.readline()    
            self.in_type_2 = str(Input.readline())
            
            Input.readline()
            self.inlet_radius_2 = float(Input.readline())
            
            Input.readline()    
            self.propellant_2 = str(Input.readline())
            
            Input.readline()    
            self.m_2 = float(Input.readline())
            
            Input.readline()
            self.delta_p_2 = float(Input.readline())
            
            Input.readline()    
            self.alpha2_2 = float(Input.readline())
            
            Input.readline()    
            self.n_2 = float(Input.readline())
            
            Input.readline()    
            self.opening_coef_2 = float(Input.readline())
            
            Input.readline()
            self.trans_angle_2 = float(Input.readline())
            
            Input.readline()    
            self.ratio_l_in_2 = float(Input.readline())
            
            Input.readline()    
            self.ratio_l_s_2 = float(Input.readline())
            
            Input.readline()    
            self.ratio_l_n_2 = float(Input.readline())
            
            Input.readline()
            self.ratio_l_n_out = float(Input.readline())
            
            Input.readline()    
            self.trans_angle_out = float(Input.readline())
            
            
            print("ST2: Data read Successful")   
            
        else:
            print("ERROR: Configuration type not recognized")
        
        
        ##### Propellant properties #####
        
        Input = open('prop_data.txt','r') 
        
        for i in range(4): Input.readline()
        
        Input.readline()
        self.n_propellants = int(Input.readline())
        
        Input.readline()
        self.n_properties = int(Input.readline())
        
        self.properties = []
        
        
        for i in range(self.n_propellants):
            
            self.properties.append([])   
            
            Input.readline()
            self.properties[i].append(str(Input.readline()))
            
            for j in range(self.n_properties-1):
        
                Input.readline()
                self.properties[i].append(float(Input.readline()))
                
            Input.readline()
        
        
        
        prop_ST1 = []
        
        prop_ST2 = []
        
        
        if(self.config == "Mono\n"):
            
            ##### properties ST1 #####
            
            for i in range(self.n_propellants):
                if(self.properties[i][0] == self.propellant_1):
                    prop_ST1 = self.properties[i]
            
            self.rho_1 = prop_ST1[1]
            
            self.din_visc_1 = prop_ST1[2]
            
            self.visc_sup_1 = prop_ST1[3]
            
        elif(self.config == "Bi\n"):
        
            ##### properties ST1 #####
            
            for i in range(self.n_propellants):
                if(self.properties[i][0] == self.propellant_1):
                    prop_ST1 = self.properties[i]
            
            self.rho_1 = prop_ST1[1]
            
            self.din_visc_1 = prop_ST1[2]
            
            self.visc_sup_1 = prop_ST1[3]
        
        
            ##### properties ST2 #####
        
            for i in range(self.n_propellants):
                if(self.properties[i][0] == self.propellant_2):
                    prop_ST2 = self.properties[i]
        
            self.rho_2 = prop_ST2[1]
            
            self.din_visc_2 = prop_ST2[2]
            
            self.visc_sup_2 = prop_ST2[3]
            
            
        else:
            print("ERROR: Configuration type not recognized (prop)")
        
