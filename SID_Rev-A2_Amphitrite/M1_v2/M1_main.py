# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 2020

@author: A. Goulart
"""

##### Data File Slection #####

#folder_name = "rev_A01_v1"
#folder_name = "rev_A01_v2"
##### ##### #####
#folder_name = "PION"
##### ##### #####
folder_name = "tests"

#file_name = "A01_FC"
#file_name = "A01_Main"
##### ##### #####
#file_name = "PION_10_F"
#file_name = "PION_10_Ox_Outer"
#file_name = "PION_10_Ox_Inner"
##### ##### #####
#file_name = "RD0110"
file_name = "RD0110_Ox"
#file_name = "Fuel_1"
#file_name = "Oxidizer_1"


from M1_Method_v2 import Method_1

Run_1 = Method_1(folder_name, file_name)

Run_1.set_it_lim(100)
Run_1.set_erro_max(1e-6)

Run_1.show_plot_error(1)
Run_1.show_plot_inj(1)

Run_1.set_dpi(200)
Run_1.set_textbox_fontsize(6)
Run_1.set_textbox_offset(24)

Run_1.run_M1()
