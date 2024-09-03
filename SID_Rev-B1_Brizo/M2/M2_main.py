# -*- coding: utf-8 -*-
"""
@author: A. Goulart
"""

##### ##### ##### ##### #####
#### Data File Selection #####
##### ##### ##### ##### #####

##### ##### ##### ##### #####
Foldername = "rev_B01/v_3"

##

#Filename = "B01_C_Ox"
#Filename = "B01_C_F"

#Filename = "B01_W_Ox"
Filename = "B01_W_F"
##### ##### ##### ##### #####


##### ##### ##### ##### ##### OK
#Foldername = "RD0110/v_2"

##

#Filename = "RD0110_Ox"
#Filename = "RD0110_F"
##### ##### ##### ##### #####


##### ##### ##### ##### #####
#Foldername = "RD0109/v_2"

##

#Filename = "RD0109_Ox"
#Filename = "RD0109_F"
##### ##### ##### ##### #####


from M2_Method import Method_2

Run_1 = Method_2(Foldername, Filename)

Run_1.set_it_lim(100)
Run_1.set_erro_max(1e-6)

Run_1.set_delta_p_i(0.50)
Run_1.set_delta_p_f(10.00)
Run_1.set_delta_p_step(0.05)

Run_1.set_lambda_type("Bazarov")
#Run_1.set_lambda_type("Bayvel")

Run_1.run_M2()
