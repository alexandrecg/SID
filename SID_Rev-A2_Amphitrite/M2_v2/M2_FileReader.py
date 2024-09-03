# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 07:40:15 2020

@author: A. Goulart
"""



class FileReader_2:
    def setFolderName(self,foname):
        self.Foldername = foname
    def setFileName(self,finame):
        self.Filename = finame
    def read(self):
        '''
        ##### Injector parameters reading #####
        
        
        
        file_loc = str(self.foldername) + '/' + str(self.filename) + '.txt'
        
        with open(str(file_loc)) as f:
            for i, l in enumerate(f):
                pass
        f_len = i + 1
        
        Input = open(str(file_loc),'r') 
    
        for i in range(5): Input.readline()
        
        Injector_ID = str(Input.readline())
        print("Identificação do injetor: \t", Injector_ID)
        
        Input.readline()
        Config = str(Input.readline())
        print("Configuração do injetor: \t", Config)
        
        print("File Length: ", f_len)
        
        DataStream = []
        
        for i in range(f_len):
            Input.readline()
            DataStream.append(Input.readline())
                
        
        
        
        
        def file_len(fname):
            with open(fname) as f:
                for i, l in enumerate(f):
                    pass
            return i + 1
        '''
        
        #from M2_main import Foldername, Filename
        
        ##### ##### ##### ##### ##### ##### #####
        #####  Injector parameters reading  #####
        ##### ##### ##### ##### ##### ##### #####
        
        print("\033[1m"+"\n Filereader initialization... \n"+"\033[0m")
        
        self.file_loc = str(self.Foldername) + '/' + str(self.Filename) + '.txt'
        
        print("\nFile location: ",self.file_loc,"\n")
        
        with open(str(self.file_loc)) as self.f:
            for self.i, self.l in enumerate(self.f):
                pass
        self.f_len = self.i + 1
        
        
        #f_len = file_len(str(file_loc))
        
        self.Input = open(str(self.file_loc),'r') 
        
        for self.i in range(5): self.Input.readline()
        
        self.Injector_ID = str(self.Input.readline())
        print("Identificação do injetor: \t", self.Injector_ID)
        
        self.Input.readline()
        self.Config = str(self.Input.readline())
        print("Configuração do injetor: \t", self.Config)
        
        print("File Length: ", self.f_len)
        
        self.DataStream = []
        
        for i in range(self.f_len):
            self.Input.readline()
            self.DataStream.append(self.Input.readline())
        
        
        if(self.Config == "Mono\n"):
            
            ##### Beginning of ST1 data #####
            
            self.t_w = float(self.DataStream[0])
            
            self.In_type_1 = str(self.DataStream[1])
            
            self.Propellant_1 = str(self.DataStream[2])
            
            self.m_1 = float(self.DataStream[3])
            
            self.Delta_P0_1 = float(self.DataStream[4])
            
            self.N_1 = int(self.DataStream[5])
            
            self.d_in_1 = float(self.DataStream[6])
    
            self.L_in_1 = float(self.DataStream[7])
            
            self.in_chamfer_1 = float(self.DataStream[8])
            
            self.d_s_1 = float(self.DataStream[9])
            
            self.L_s_1 = float(self.DataStream[10])
            
            self.d_n_1 = float(self.DataStream[11])
            
            self.L_n_1 = float(self.DataStream[12])
            
            self.R_0_1 = float(self.DataStream[13])
            
        elif(self.Config == "Bi\n"):
        
            ##### Beginning of ST1 data #####
        
            self.t_w = float(self.DataStream[0])
            
            self.In_type_1 = str(self.DataStream[1])
            
            self.Propellant_1 = str(self.DataStream[2])
            
            self.m_1 = float(self.DataStream[3])
            
            self.Delta_P0_1 = float(self.DataStream[4])
            
            self.N_1 = int(self.DataStream[5])
            
            self.d_in_1 = float(self.DataStream[6])
    
            self.L_in_1 = float(self.DataStream[7])
            
            self.in_chamfer_1 = float(self.DataStream[8])
            
            self.d_s_1 = float(self.DataStream[9])
            
            self.L_s_1 = float(self.DataStream[10])
            
            self.d_n_1 = float(self.DataStream[11])
            
            self.L_n_1 = float(self.DataStream[12])
            
            self.R_0_1 = float(self.DataStream[13])
            
            
            ##### Beginning of ST2 data #####
            
            self.In_type_2 = str(self.DataStream[14])
            
            self.Propellant_2 = str(self.DataStream[15])
            
            self.m_2 = float(self.DataStream[16])
            
            self.Delta_P0_2 = float(self.DataStream[17])
            
            self.N_2 = int(self.DataStream[18])
            
            self.d_in_2 = float(self.DataStream[19])
            
            self.L_in_2 = float(self.DataStream[20])
            
            self.in_chamfer_2 = float(self.DataStream[21])
            
            self.d_s_2 = float(self.DataStream[22])
            
            self.L_s_2 = float(self.DataStream[23])
            
            self.d_n_2 = float(self.DataStream[24])
            
            self.L_n_2 = float(self.DataStream[25])
            
            self.R_0_2 = float(self.DataStream[26])
            
            self.R_out_2 = float(self.DataStream[27])
            
            
        else:
            print("ERROR: Configuration type not recognized")
        
        
        
        
        
        ##### ##### ##### ##### ##### #####
        #####  Propellant properties  #####
        ##### ##### ##### ##### ##### #####
        
        self.Input = open('prop_data.txt','r') 
        
        for i in range(4): self.Input.readline()
        
        self.Input.readline()
        self.n_propellants = int(self.Input.readline())
        
        self.Input.readline()
        self.n_properties = int(self.Input.readline())
        
        self.Properties = []
        
        
        for i in range(self.n_propellants):
            
            self.Properties.append([])   
            
            self.Input.readline()
            self.Properties[i].append(str(self.Input.readline()))
            
            for j in range(self.n_properties-1):
        
                self.Input.readline()
                self.Properties[i].append(float(self.Input.readline()))
                
            self.Input.readline()
        
        
        
        self.prop_ST1 = []
        
        self.prop_ST2 = []
        
        
        if(self.Config == "Mono\n"):
            
            ##### properties ST1 #####
            
            for i in range(self.n_propellants):
                if(self.Properties[i][0] == self.Propellant_1):
                    self.prop_ST1 = self.Properties[i]
            
            self.Rho_1 = self.prop_ST1[1]
            
            self.Din_visc_1 = self.prop_ST1[2]
            
            self.Visc_sup_1 = self.prop_ST1[3]
            
        elif(self.Config == "Bi\n"):
        
            ##### properties ST1 #####
            
            for i in range(self.n_propellants):
                if(self.Properties[i][0] == self.Propellant_1):
                    self.prop_ST1 = self.Properties[i]
            
            self.Rho_1 = self.prop_ST1[1]
            
            self.Din_visc_1 = self.prop_ST1[2]
            
            self.Visc_sup_1 = self.prop_ST1[3]
        
        
            ##### properties ST2 #####
        
            for i in range(self.n_propellants):
                if(self.Properties[i][0] == self.Propellant_2):
                    self.prop_ST2 = self.Properties[i]
        
            self.Rho_2 = self.prop_ST2[1]
            
            self.Din_visc_2 = self.prop_ST2[2]
            
            self.Visc_sup_2 = self.prop_ST2[3]
            
            
        else:
            print("ERROR: Configuration type not recognized (prop)")
        
            

        print("\033[1m"+"\n File reading complete!  \n\n"+"\033[0m")