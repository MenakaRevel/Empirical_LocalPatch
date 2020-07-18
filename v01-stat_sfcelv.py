#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import sys
import os
import errno
from numpy import ma 
import matplotlib.gridspec as gridspec
import string
from scipy.fftpack import fft, ifft, fftfreq
#from slacker import Slacker
from multiprocessing import Pool
from multiprocessing import Process
import xarray as xr

import params as pm
#--
#----
def mk_dir(sdir):
    try:
        os.makedirs(sdir)
    except:
        pass
#--
mk_dir(pm.out_dir()+"/stats")
#--read outflow netCDF4--
tag="%04d-%04d"%(pm.starttime()[0],pm.endtime()[0])
# sfcelv
fname=pm.out_dir()+"/CaMa_out/"+pm.input_name()+"/sfcelv"+tag+".nc"
nc=xr.open_dataset(fname)
#
#sfcelv_mean=nc.sfcelv.mean(axis=0)
#sfcelv_std=nc.sfcelv.std(axis=0)
#sfcelv_mean=sfcelv_mean.values
#sfcelv_std=sfcelv_std.values
#sfcelv_mean.tofile(pm.out_dir()+"/stats/sfcelv_mean"+tag+".bin")
#sfcelv_std.tofile(pm.out_dir()+"/stats/sfcelv_std"+tag+".bin")
#print (np.shape(sfcelv_mean))
#-----
year=2003
tag="%04d-%04d"%(year,year)
#data=nc.sfcelv.loc('2003-01-01','2003-12-31','sfcelv')
syear=pm.starttime()[0]
smon=pm.starttime()[1]
sday=pm.starttime()[2]
st=datetime.date(syear,smon,sday)
d1=datetime.date(year,1,1)
d2=datetime.date(year,12,31)
print st
dt1=(d1-st).days - 1
dt2=(d2-st).days
print dt1, dt2
data=nc.sfcelv[dt1:dt2]
print data, np.shape(data.values)
sfcelv_mean=data.mean(axis=0)
sfcelv_std=data.std(axis=0)
sfcelv_mean=sfcelv_mean.values
sfcelv_std=sfcelv_std.values
print (np.shape(sfcelv_mean))
sfcelv_mean.tofile(pm.out_dir()+"/stats/sfcelv_mean"+tag+".bin")
sfcelv_std.tofile(pm.out_dir()+"/stats/sfcelv_std"+tag+".bin")
nc.close()
