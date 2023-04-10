#!/opt/local/bin/python
# -*- coding: utf-8 -*-
""" make plots for fitted semivariograms
for each upstream and downstream
Menaka@IIS 2023/04/06"""
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import matplotlib as mpl
import sys
import os
import matplotlib.gridspec as gridspec
import string
import calendar
import errno
import re
import math
from numpy import ma 
import matplotlib.gridspec as gridspec
#from mpl_toolkits.basemap import Basemap
#from slacker import Slacker
from multiprocessing import Pool
from multiprocessing import Process
#---
import LMA_semivari as LMA
#---
# fname="../semivar/conus_06min_VIC_BC/03910227/up00312.svg"
fname="../semivar/glb_15min_S14FD/03570231/up00046.svg"
with open(fname, "r") as f:
    lines = f.readlines()
ldis=[]
lgamma=[]
for line in lines[1::]:
    print (line)
    line = filter(None, re.split(" ",line))
    dis  = float(line[2])
    gamma= float(line[3])
    ldis.append(dis)
    lgamma.append(gamma)
# plot
fig = plt.figure()
ax  = fig.add_subplot(111)
ax.plot(ldis,lgamma,marker="o",linestyle="None",markeredgecolor="k",markerfacecolor=None,markersize=5)
plt.show()