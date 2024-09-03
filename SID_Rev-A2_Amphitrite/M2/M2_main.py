# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 07:40:55 2020

@author: A. Goulart
"""

##### Data File Slection #####

Foldername = "rev_A01"
#Foldername = "REF"

#Filename = "A01_FC_orig"
#Filename = "A01_FC_mod"
#Filename = "A01_Main_orig"
Filename = "A01_Main_mod"
#Filename = "RD0110"
#Filename = "ALVES"

from M2_Method import Method_2

Run_1 = Method_2(Foldername, Filename)

Run_1.run_M2()
