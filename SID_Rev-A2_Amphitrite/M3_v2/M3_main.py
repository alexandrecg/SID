# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 07:40:55 2020

@author: A. Goulart
"""

##### Data File Slection #####


Foldername = "tests"


Filename = "RD0110_Ox"


from M3_Method_v2 import Method_3

Run_1 = Method_3(Foldername, Filename)

Run_1.set_it_lim(100)
Run_1.set_erro_max(1e-6)
Run_1.set_sample_size(10)

Run_1.set_delta_p_i(0.50)
Run_1.set_delta_p_f(10.00)
Run_1.set_delta_p_step(0.05)

Run_1.run_M3()
