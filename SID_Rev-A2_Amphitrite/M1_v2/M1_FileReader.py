# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 2020

@author: A. Goulart
"""

class FileReader_1:
    def setFolderName(self,foname):
        self.foldername = foname
    def setFileName(self,finame):
        self.filename = finame
    def read(self):

        ##### Injector parameters reading #####
        
        file_loc = str(self.foldername) + '/' + str(self.filename) + '.txt'
        
        Input = open(str(file_loc),'r') 
        
        for i in range(5): Input.readline()
        
        Injector_ID = str(Input.readline())
        print("Identificação do injetor: \t", Injector_ID)
        
        Input.readline()
        self.Config = str(Input.readline())
        print("Configuração do injetor: \t", self.Config)
        
        
        if(self.Config == "Mono\n"):
            
            ##### Beginning of ST1 data #####
            
            Input.readline()
            self.t_w = float(Input.readline())
            
            Input.readline()
            self.delta = float(Input.readline())
            
            Input.readline()
            self.angle_dif = float(Input.readline())
            
            Input.readline()
            self.recess = float(Input.readline())
            
            
            ##### Beginning of ST1 data #####
            
            Input.readline()    
            self.In_type_1 = str(Input.readline())
            
            Input.readline()
            self.entrance_radius_1 = float(Input.readline())
            
            Input.readline()    
            self.Propellant_1 = str(Input.readline())
            
            Input.readline()    
            self.m_1 = float(Input.readline())
            
            Input.readline()    
            self.Delta_P_1 = float(Input.readline())
            
            Input.readline()    
            self.Alpha2_1 = float(Input.readline())
            
            Input.readline()    
            self.n_1 = int(Input.readline())
            
            Input.readline()    
            self.Rs_Rn_1 = float(Input.readline())
            
            Input.readline()
            self.Trans_angle_1 = float(Input.readline())
            
            Input.readline()    
            self.LD_in_1 = float(Input.readline())
            
            Input.readline()    
            self.LD_s_1 = float(Input.readline())
            
            Input.readline()    
            self.LD_n_1 = float(Input.readline())
            
            print("Data read Successful")        
            
        elif(self.Config == "Bi\n"):
            
            Input.readline()
            self.t_w = float(Input.readline())
            
            Input.readline()
            self.delta = float(Input.readline())
            
            Input.readline()
            self.angle_dif = float(Input.readline())
            
            Input.readline()
            self.recess = float(Input.readline())
            
            
            ##### Beginning of ST1 data #####
            
            Input.readline()    
            self.In_type_1 = str(Input.readline())
            
            Input.readline()
            self.entrance_radius_1 = float(Input.readline())
            
            Input.readline()    
            self.Propellant_1 = str(Input.readline())
            
            Input.readline()    
            self.m_1 = float(Input.readline())
            
            Input.readline()    
            self.Delta_P_1 = float(Input.readline())
            
            Input.readline()    
            self.Alpha2_1 = float(Input.readline())
            
            Input.readline()    
            self.n_1 = int(Input.readline())
            
            Input.readline()    
            self.Rs_Rn_1 = float(Input.readline())
            
            Input.readline()
            self.Trans_angle_1 = float(Input.readline())
            
            Input.readline()    
            self.LD_in_1 = float(Input.readline())
            
            Input.readline()    
            self.LD_s_1 = float(Input.readline())
            
            Input.readline()    
            self.LD_n_1 = float(Input.readline())
            
            ##### Beginning of ST2 data #####
            
            Input.readline()    
            self.In_type_2 = str(Input.readline())
            
            Input.readline()
            self.entrance_radius_2 = float(Input.readline())
            
            Input.readline()    
            self.Propellant_2 = str(Input.readline())
            
            Input.readline()    
            self.m_2 = float(Input.readline())
            
            Input.readline()
            self.Delta_P_2 = float(Input.readline())
            
            Input.readline()    
            self.Alpha2_2 = float(Input.readline())
            
            Input.readline()    
            self.n_2 = float(Input.readline())
            
            Input.readline()    
            self.Rs_Rn_2 = float(Input.readline())
            
            Input.readline()
            self.Trans_angle_2 = float(Input.readline())
            
            Input.readline()    
            self.LD_in_2 = float(Input.readline())
            
            Input.readline()    
            self.LD_s_2 = float(Input.readline())
            
            Input.readline()    
            self.LD_n_2 = float(Input.readline())
            
            Input.readline()
            self.LD_n_3 = float(Input.readline())
            
            Input.readline()    
            self.Trans_angle_3 = float(Input.readline())
            
            
            print("Data read Successful")   
            
        else:
            print("ERROR: Configuration type not recognized")
        
        
        ##### Propellant properties #####
        
        Input = open('prop_data.txt','r') 
        
        for i in range(4): Input.readline()
        
        Input.readline()
        self.n_propellants = int(Input.readline())
        
        Input.readline()
        self.n_properties = int(Input.readline())
        
        self.Properties = []
        
        
        for i in range(self.n_propellants):
            
            self.Properties.append([])   
            
            Input.readline()
            self.Properties[i].append(str(Input.readline()))
            
            for j in range(self.n_properties-1):
        
                Input.readline()
                self.Properties[i].append(float(Input.readline()))
                
            Input.readline()
        
        
        
        prop_ST1 = []
        
        prop_ST2 = []
        
        
        if(self.Config == "Mono\n"):
            
            ##### properties ST1 #####
            
            for i in range(self.n_propellants):
                if(self.Properties[i][0] == self.Propellant_1):
                    prop_ST1 = self.Properties[i]
            
            self.Rho_1 = prop_ST1[1]
            
            self.Din_visc_1 = prop_ST1[2]
            
            self.Visc_sup_1 = prop_ST1[3]
            
        elif(self.Config == "Bi\n"):
        
            ##### properties ST1 #####
            
            for i in range(self.n_propellants):
                if(self.Properties[i][0] == self.Propellant_1):
                    prop_ST1 = self.Properties[i]
            
            self.Rho_1 = prop_ST1[1]
            
            self.Din_visc_1 = prop_ST1[2]
            
            self.Visc_sup_1 = prop_ST1[3]
        
        
            ##### properties ST2 #####
        
            for i in range(self.n_propellants):
                if(self.Properties[i][0] == self.Propellant_2):
                    prop_ST2 = self.Properties[i]
        
            self.Rho_2 = prop_ST2[1]
            
            self.Din_visc_2 = prop_ST2[2]
            
            self.Visc_sup_2 = prop_ST2[3]
            
            
        else:
            print("ERROR: Configuration type not recognized (prop)")
        
