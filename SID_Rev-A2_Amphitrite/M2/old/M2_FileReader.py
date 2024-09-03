# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 07:40:15 2020

@author: A. Goulart
"""
'''
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
'''

from M2_main import Foldername, Filename

##### ##### ##### ##### ##### ##### #####
#####  Injector parameters reading  #####
##### ##### ##### ##### ##### ##### #####

file_loc = str(Foldername) + '/' + str(Filename) + '.txt'


with open(str(file_loc)) as f:
    for i, l in enumerate(f):
        pass
f_len = i + 1


#f_len = file_len(str(file_loc))

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

print(DataStream)

if(Config == "Mono\n"):
    
    ##### Beginning of ST1 data #####
    '''   
    N_in = 3

    print("Número de canais de entrada (N_in): \t\t%.1f"%(N_in))
    
    d_in = 0.80 #mm
    
    print("Diâmetro dos canais de entrada (d_in): \t\t%.2f \tmm"%(d_in))
    
    r_in = d_in/2 #mm
    
    k_L = 2.0 #2.0 a 3.0
    
    #d_bx = 1.1* d_in#0.20 #mm (Chamfro de entrada)
    
    L_in = k_L*d_in
    
    print("Comprimento dos canais de entrada (L_in): \t%.2f \tmm"%(L_in))
    
    D_cv = 4.00 #mm
    
    print("Diâmetro da Câmara de Vórtice (D_cv): \t\t%.2f \tmm"%(D_cv))
    
    R_cv = D_cv/2 #mm
    
    L_cv = 6 #mm
    
    print("Comprimento da Câmara de Vórtice (L_cv): \t%.2f \tmm"%(L_cv))
    
    D_n = 2.00 #mm ###################################
    
    print("Diâmetro do Orifício de saída (D_n): \t\t%.2f \tmm"%(D_n))
    
    R_n = D_n/2 #mm
    
    R_in = R_cv-r_in-0.05 #mm ##################################
    
    print("Raio de Entrada dos canais (R_in): \t\t%.2f \tmm"%(R_in))
    
    A_n = np.pi*R_n**2
    
    

    
    
    '''
   
elif(Config == "Bi\n"):
    '''
    Input.readline()
    t_w = float(Input.readline())
    
    Input.readline()
    delta = float(Input.readline())
    
    Input.readline()
    angle_dif = float(Input.readline())
    
    Input.readline()
    recess = float(Input.readline())
    '''
    
    ##### Beginning of ST1 data #####
    
    
    
    ##### Beginning of ST2 data #####
    
    

else:
    print("ERROR: Configuration type not recognized")





##### ##### ##### ##### ##### #####
#####  Propellant properties  #####
##### ##### ##### ##### ##### #####

Input = open('prop_data.txt','r') 

for i in range(4): Input.readline()

Input.readline()
n_propellants = int(Input.readline())

Input.readline()
n_properties = int(Input.readline())

Properties = []


for i in range(n_propellants):
    
    Properties.append([])   
    
    Input.readline()
    Properties[i].append(str(Input.readline()))
    
    for j in range(n_properties-1):

        Input.readline()
        Properties[i].append(float(Input.readline()))
        
    Input.readline()



prop_ST1 = []

prop_ST2 = []


if(Config == "Mono\n"):
    
    ##### properties ST1 #####
    
    for i in range(n_propellants):
        if(Properties[i][0] == Propellant_1):
            prop_ST1 = Properties[i]
    
    Rho_1 = prop_ST1[1]
    
    Din_visc_1 = prop_ST1[2]
    
    Visc_sup_1 = prop_ST1[3]
    
elif(Config == "Bi\n"):

    ##### properties ST1 #####
    
    for i in range(n_propellants):
        if(Properties[i][0] == Propellant_1):
            prop_ST1 = Properties[i]
    
    Rho_1 = prop_ST1[1]
    
    Din_visc_1 = prop_ST1[2]
    
    Visc_sup_1 = prop_ST1[3]


    ##### properties ST2 #####

    for i in range(n_propellants):
        if(Properties[i][0] == Propellant_2):
            prop_ST2 = Properties[i]

    Rho_2 = prop_ST2[1]
    
    Din_visc_2 = prop_ST2[2]
    
    Visc_sup_2 = prop_ST2[3]
    
    
else:
    print("ERROR: Configuration type not recognized (prop)")

    
    
    