# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:33:36 2020

@author: A. Goulart
"""

import matplotlib.pyplot as plt


fig, ax = plt.subplots(figsize=(3, 3))

ax.annotate("Test",xy=(0.2, 0.2), xycoords='data',xytext=(0.8, 0.8))

plt.show()