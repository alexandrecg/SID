# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 2020

@author: A. Goulart
"""

##### Data File Slection #####

Foldername = "rev_A01"
##### ##### #####
#Foldername = "PION"

Filename = "A01_FC"
#Filename = "A01_Main"
##### ##### #####
#Filename = "PION_10_F"
#Filename = "PION_10_Ox_Outer"
#Filename = "PION_10_Ox_Inner"



from M1_Method_v15 import Method_1

Run_1 = Method_1(Foldername, Filename)

#Run_1.set_it_lim(10000)
Run_1.set_erro_max(1e-6)

Run_1.run_M1()
