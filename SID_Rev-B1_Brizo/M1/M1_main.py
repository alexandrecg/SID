# -*- coding: utf-8 -*-
"""
@author: A. Goulart
"""

##### ##### ##### ##### #####
folder_name = "rev_B01/v_2"
##### ##### #####
#file_name = "B01_C"
file_name = "B01_W"
##### ##### ##### ##### #####


##### ##### ##### ##### #####
#folder_name = "RD0110/v_2"
##### ##### #####
#file_name = "RD0110"
##### ##### ##### ##### #x####

##### ##### ##### ##### #####
#folder_name = "RD0109"
##### ##### #####
#file_name = "RD0109"
##### ##### ##### ##### #####


from M1_Method import Method_1

Run_1 = Method_1(folder_name, file_name)

Run_1.set_it_lim(100)
Run_1.set_erro_max(1e-6)

Run_1.set_show_it(0)
Run_1.show_plot_error(0)
Run_1.show_plot_inj(1)

Run_1.set_dpi(200)
Run_1.set_textbox_fontsize(6)
#Run_1.set_textbox_offset(24)
#Run_1.set_textbox_offset(6)
Run_1.set_textbox_offset(14)

Run_1.run_M1()
