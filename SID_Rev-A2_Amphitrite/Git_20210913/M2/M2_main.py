# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 07:40:55 2020

@author: A. Goulart
"""

##### Data File Slection #####

#Foldername = "rev_A01"
##### ##### #####
Foldername = "tests"


#Filename = "A01_FC_orig"
#Filename = "A01_FC_mod"
#Filename = "A01_Main_orig"
#Filename = "A01_Main_mod"
##### ##### #####
#Filename = "RD0110_Ox"
#Filename = "RD0110_Ox_M1"
Filename = "RD0110_Ox_mod"


from M2_Method_v2 import Method_2

Run_1 = Method_2(Foldername, Filename)

Run_1.set_it_lim(100)
Run_1.set_erro_max(1e-6)

Run_1.set_delta_p_i(0.50)
Run_1.set_delta_p_f(10.00)
Run_1.set_delta_p_step(0.05)

Run_1.run_M2()
