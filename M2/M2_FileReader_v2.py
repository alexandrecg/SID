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
         
        ##### ##### ##### ##### ##### ##### #####
        #####  Injector parameters reading  #####
        ##### ##### ##### ##### ##### ##### #####
        
        print("\033[1m"+"\n Filereader initialization... \n"+"\033[0m")
        
        self.file_loc = str(self.Foldername) + '/' + str(self.Filename) + '.txt'
        
        print("\nFile location: ",self.file_loc,"\n")

 
        self.f_len = sum(1 for line in open(str(self.file_loc)))
            

        self.Input = open(str(self.file_loc),'r') 
        
        for i in range(5): self.Input.readline()
        
        self.Injector_ID = str(self.Input.readline())
        print("Identificação do injetor: \t", self.Injector_ID)
        
        print("File Length: ", self.f_len)
        
        self.DataStream = []
        
        for i in range(self.f_len):
            self.Input.readline()
            self.DataStream.append(self.Input.readline())
        
        
            
        ##### Beginning of ST1 data #####
        
        self.inlet_type = str(self.DataStream[0])
        
        self.propellant = str(self.DataStream[1])
        
        self.m_0 = float(self.DataStream[2])
        
        self.delta_p0 = float(self.DataStream[3])
        
        self.n = int(self.DataStream[4])
        
        self.r_in_orf = float(self.DataStream[5])

        self.l_in = float(self.DataStream[6])
        
        self.inlet_radius = float(self.DataStream[7])
        
        self.r_s = float(self.DataStream[8])
        
        self.l_s = float(self.DataStream[9])
        
        self.r_n = float(self.DataStream[10])
        
        self.l_n = float(self.DataStream[11])
        
        self.r_in_pos = float(self.DataStream[12])
            
        
        
        
        
        
        
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
            
        ##### properties ST1 #####
        
        for i in range(self.n_propellants):
            if(self.Properties[i][0] == self.propellant):
                self.prop_ST1 = self.Properties[i]
        
        self.rho = self.prop_ST1[1]
        
        self.din_visc = self.prop_ST1[2]
        
        self.visc_sup = self.prop_ST1[3]
            
        
            

        print("\033[1m"+"\n File reading complete!  \n\n"+"\033[0m")