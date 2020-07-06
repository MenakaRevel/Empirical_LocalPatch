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

argvs = sys.argv

#----
def SWOT_day(yyyy,mm,dd):
  st_year,st_month,st_date=pm.starttime()
  start_time=datetime.date(st_year,st_month,st_date)
  this_time=datetime.date(int(yyyy),int(mm),int(dd))
  days=this_time-start_time
  days=days.days
  return days%21+1
#---
def slink(src,dst):
  try:
    os.symlink(src,dst)
  except OSError as exc:  # Python >2.5
    if exc.errno == errno.EEXIST:
      os.remove(dst)
      os.symlink(src,dst)
    else:
      raise
#----
#----
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#--
def fft_90(y,coeff,T,n):
# FFT of given timeseries
# Inverse FFT of (1-coeff)% of the timeseries
# get a timeseries only frequencies which have streagth geater than (1-coeff)%
# T is the sampling perios xf gives in cycles per T
  # fft of y
  yf=fft(y)
  # frequncies corresponds to yf
  xf=fftfreq(n,T)
  # maximum fft stregth
  yfmax=np.max(np.abs(yf[1::]))
  yf1=yf.copy() 
  # remove fft stregths lower than coeff*yfmax 
  yf1[1::]=ma.masked_where(np.abs(yf[1::])<coeff*yfmax,yf1[1::]).filled(0.0)
  # inverse fft of frequencies greater than (1-coeff)*yfmax 
  yif =ifft(yf1)   
  return yif,yf,xf
#--
def detrend(x,y,n):
# assume trend line as y = a + bx
  x = np.array(x)
  y = np.array(y)
  # Exy
  Exy = sum(x*y)
  Ex  = sum(x)
  Ey  = sum(y)
  Ex2 = sum(x**2)
  #--
  b = ((n*Exy) - (Ex*Ey))/((n*Ex2) - Ex**2)
  a = np.mean(y) - b*np.mean(x)
  return a + b*(x)  
#---
#def sendslack(msg):
#  slack=Slacker('xoxb-316041665447-hNXDYiTz9PZMITg0gfRGzyrC')
#  message="` "+os.path.basename(__file__)+"` : "+msg
#  slack.chat.post_message('#camasim',message);
#  return 0
#---     
#--
#slink("/hydro/covariance/CaMa_simulation/params.py", "params.py")
#slink("/hydro/SWOTDA/img_code/read_grdc.py","read_grdc.py")
import params as pm
import read_grdc as grdc
#--
mk_dir(pm.out_dir()+"/figures")
mk_dir(pm.out_dir()+"/figures/timeseries")
#--read outflow netCDF4--
tag="%04d-%04d"%(pm.starttime()[0],pm.endtime()[0])
# sfcelv
fname=pm.out_dir()+"/CaMa_out/"+pm.input_name()+"/sfcelv"+tag+".nc"
nc_sfcelv=xr.open_dataset(fname)
# removed trend line
fname=pm.out_dir()+"/CaMa_out/"+pm.input_name()+"/rmdtrnd"+tag+".nc"
nc_rmdtrend=xr.open_dataset(fname)
# removed seasonality
fname=pm.out_dir()+"/CaMa_out/"+pm.input_name()+"/rmdsesn"+tag+".nc"
nc_rmdseson=xr.open_dataset(fname)
# standardized
fname=pm.out_dir()+"/CaMa_out/"+pm.input_name()+"/standardized"+tag+".nc"
nc_standz=xr.open_dataset(fname)

#--
syear=pm.starttime()[0]
smonth=pm.starttime()[1]
sdate=pm.starttime()[2]
start_dt=datetime.date(syear,smonth,sdate)
size=60

start=0
#last_dt=datetime.date(int(argvs[1]),int(argvs[2]),int(argvs[3]))
#eyear,emonth,edate=pm.patch_end()
eyear=pm.endtime()[0]
emonth=pm.endtime()[1]
edate=pm.endtime()[2]

last_dt=datetime.date(eyear,emonth,edate)
last=(last_dt-start_dt).days + 1
N=int(last)
green2="greenyellow"

#sendslack(" _starts_ ")
#colors = pc.color_assimd()
pname=[]
xlist=[]
ylist=[]
river=[]
staid=[]
#--
#rivernames  = ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS","INDUS"]
rivernames = ["AMAZON"]#"MEKONG"]
#rivernames = grdc.grdc_river_name()
for rivername in rivernames:
  path = pm.out_dir()+"/figures/timeseries/%s"%(rivername)
  print path
  mk_dir(path)
  grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername)
  print rivername, station_loc
  river.append([rivername]*len(station_loc))
  pname.append(station_loc)
  staid.append(grdc_id)
  xlist.append(x_list)
  ylist.append(y_list)
#--
#for rivername in ["LENA","NIGER","INDUS","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS"]:
#for rivername in ["MEKONG","AMAZONAS"]:
###  if rivername=="LENA":
###      pname.append(["L1","L2","L3"])
###      xlist.append([1233,1218,1206])
###      ylist.append([  71,  98, 117])
###      river.append([rivername,rivername,rivername])
###  if rivername=="NIGER":
###      pname.append(["N1","N2","N3","N4","N5","N6"])
###      xlist.append([744,744,732,712,704,700])
###      ylist.append([342,324,310,292,295,303])
###      river.append([rivername,rivername,rivername,rivername,rivername,rivername])
###  if rivername=="AMAZONAS":
###      pname.append(["B","E","F","G"])
###      xlist.append([515,447,464,420])
###      ylist.append([364,366,416,367])
###      river.append([rivername,rivername,rivername,rivername])
###  if rivername=="MEKONG":
###      pname.append(["Me_D","Me_M","Me_U"])
###      xlist.append([1143,1139,1127])
###      ylist.append([ 319, 293, 282])
###      river.append([rivername,rivername,rivername])
###  if rivername=="MISSISSIPPI":
###      pname.append(["Mi_D","Mi_M","Mi_U","Mi_U2"]) #最後はミズーリ川
###      xlist.append([361,362,345,308])
###      ylist.append([241,214,137,167])
###      river.append([rivername,rivername,rivername,rivername])
###  if rivername=="OB":
###      pname.append(["O_D","O_M","O_U"])
###      xlist.append([995,996,1048])
###      ylist.append([ 92,121, 159])
###      river.append([rivername,rivername,rivername])
###  if rivername=="CONGO":
###      pname.append(["C_D","C_M","C_U"])
###      xlist.append([772,813,834])
###      ylist.append([383,353,397])
###      river.append([rivername,rivername,rivername])
###  if rivername=="INDUS":
###      pname.append(["I_D","I_Sub","I_M"])
###      xlist.append([992,995,1003])
###      ylist.append([251,252,233])
###      river.append([rivername,rivername,rivername])

river=([flatten for inner in river for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])


pnum=len(pname)

org=[]
tdl=[]
fft=[]
sdz=[]
#swt={}
#for point in np.arange(pnum):
#    swt[point] = []

###for day in np.arange(start,last):
###    target_dt=start_dt+datetime.timedelta(days=day)
###    yyyy='%04d' % (target_dt.year)
###    mm='%02d' % (target_dt.month)
###    dd='%02d' % (target_dt.day)
###    print yyyy,mm,dd
###
####    if target_dt.month==1 and target_dt.day==1:
####      sendslack(" _reading_ "+yyyy)
###
###
###    # make org
###    fname=pm.out_dir()+"/CaMa_out/"+yyyy+"/"+yyyy+mm+dd+"/sfcelv"+yyyy+".bin"
###    orgfile=np.fromfile(fname,np.float32).reshape([720,1440])
###
###    # make tdl
###    fname=pm.out_dir()+"/CaMa_out/"+yyyy+"/"+yyyy+mm+dd+"/tdline"+yyyy+".bin"
###    tdlfile=np.fromfile(fname,np.float32).reshape([720,1440])
###
###    # make fft
###    fname=pm.out_dir()+"/CaMa_out/"+yyyy+"/"+yyyy+mm+dd+"/sesanlty"+yyyy+".bin"
###    fftfile=np.fromfile(fname,np.float32).reshape([720,1440])
###
###    # make standadize
###    fname=pm.out_dir()+"/CaMa_out/"+yyyy+"/"+yyyy+mm+dd+"/standze"+yyyy+".bin"
###    sdzfile=np.fromfile(fname,np.float32).reshape([720,1440])
###
###    org_frag=[]
###    tdl_frag=[]
###    fft_frag=[]
###    sdz_frag=[]
###    for point in np.arange(pnum):
###        xpoint=xlist[point]
###        ypoint=ylist[point]
###        org_frag.append(orgfile[ypoint,xpoint])
###        tdl_frag.append(tdlfile[ypoint,xpoint])
###        fft_frag.append(fftfile[ypoint,xpoint]) 
###        sdz_frag.append(sdzfile[ypoint,xpoint])
###    org.append(org_frag)
###    tdl.append(tdl_frag)
###    fft.append(fft_frag) 
###    sdz.append(sdz_frag)
###
###org=np.array(org)
###tdl=np.array(tdl)
###fft=np.array(fft)
###sdz=np.array(sdz)
####print np.shape(org)
hgt=12 
wdt=12
days=np.arange(start_dt.year,last_dt.year+1,5)
xxlist=np.linspace(1,N,len(days))
#for point in np.arange(pnum):
def make_fig(point):
    plt.close()

    #print np.shape(org[:,point])
    #print np.shape(fft[:,point]) 
    fig=plt.figure(figsize=(wdt,hgt))
    #fig=plt.figure()
    G = gridspec.GridSpec(4,1)
    # standardized data
    sdz=nc_standz.standardize[:,ylist[point],xlist[point]]
    # actuall curve
    org=nc_sfcelv.sfcelv[:,ylist[point],xlist[point]]
    # seasonality removed
    fft=nc_rmdseson.rmdsesn[:,ylist[point],xlist[point]]
    # trend line
    tdl=nc_rmdtrend.rmdtrnd[:,ylist[point],xlist[point]]
    # original data
    ax1 = fig.add_subplot(G[3,0])
    ax1.plot(np.arange(start,last),sdz,label="standardize",color="black",linestyle="-",linewidth=0.7,zorder=101)
    #ax1.plot(np.arange(start,last),org-fft,label="fft",color="green",linewidth=0.7,zorder=101)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('standardized\n data', color='k')
    ax1.set_xlabel('$year$', color='k')
    #ax1.get_yaxis().set_ticklabels([])
    ax1.set_xlim(xmin=0,xmax=last+1)
    ax1.tick_params('y', colors='k')
    ax1.set_xticks(xxlist)
    ax1.set_xticklabels(days,fontsize=8) 
    #--trend line
    ax2 = fig.add_subplot(G[1,0])
    ax2.plot(np.arange(start,last),tdl,label="trend line",color="grey",linewidth=0.7,zorder=103)
    # sesonality
    #ax2.plot(np.arange(start,last),yif,label="fft",color="r",linewidth=0.7,zorder=102)
    ax2.plot(np.arange(start,last),tdl-fft,label="seasonality",color="red",linewidth=0.7,zorder=101)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax2.set_ylabel('WSE ($m$)', color='k')
    ax2.set_xlabel('$year$)', color='k')
    #ax2.get_yaxis().set_ticklabels([])
    ax2.set_xlim(xmin=0,xmax=last+1)
    ax2.tick_params('y' , colors='k')
    ax2.set_xticks(xxlist)
    ax2.set_xticklabels(days,fontsize=8) 
    #########################
    #--senality removed data 
    ax3 = fig.add_subplot(G[2,0])
    ax3.plot(np.arange(start,last),fft,label="seasonality removed",color="xkcd:ocean",linewidth=0.7,zorder=102)
#    ax3.scatter(np.arange(start,last),fft[:,point],marker="*",color="r")
    # actual curve
    #ax3.plot(np.arange(start,last),org-tdl,label="true",color="gray",linewidth=0.7,zorder=101)
    #ax3.plot(xf[1:N//2],1.0/N * np.abs(yf[0:N//2])[1:],color="blue",linewidth=0.7,zorder=101)
    #ax3.axhline((1.0/N) * yfmax * coeff,color="g",linestyle="--",linewidth=0.7,zorder=100)
    #ax3.plot(xf[1:],psd,color="r",linewidth=0.7,zorder=101)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax3.set_ylabel('$m$', color='k')
    ax3.set_xlabel('$year$', color='k')
    ax3.set_xlim(xmin=0,xmax=last+1)
    #ax3.set_ylim(ymin=0)
    #ax2.tick_params('y', colors='k')
    ax3.set_xticks(xxlist)
    ax3.set_xticklabels(days,fontsize=8) 
    #--CaMa-Flood WSE 
    ax4 = fig.add_subplot(G[0,0])
    # actual curve
    ax4.plot(np.arange(start,last),org,label="WSE",color="b",linewidth=0.7,zorder=102)
#    ax3.scatter(np.arange(start,last),fft[:,point],marker="*",color="r")
    ax4.plot(np.arange(start,last),org-tdl,label="trend line",color="g",linewidth=0.7,zorder=101)
    # actual curve
    #ax4.plot(np.arange(start,last),org-tdl,label="trend line",color="gray",linewidth=0.7,zorder=101)
    #ax4.plot(xf[1:N//2],1.0/N * np.abs(yf[0:N//2])[1:],color="blue",linewidth=0.7,zorder=101)
    #ax4.axhline((1.0/N) * yfmax * coeff,color="g",linestyle="--",linewidth=0.7,zorder=100)
    #ax4.plot(xf[1:],psd,color="r",linewidth=0.7,zorder=101)
    # Make the y-axis label, ticks and tick labels match the line color.
    ax4.set_ylabel('$WSE$ ($m$)', color='k')
    ax4.set_xlabel('$year$', color='k')
    ax4.set_xlim(xmin=0,xmax=last+1)
    #ax4.set_ylim(ymin=0)
    #ax4.tick_params('y', colors='k')
    ax4.set_xticks(xxlist)
    ax4.set_xticklabels(days,fontsize=8) 
    #------------------
    #xff1.tofile("../img/fft_sfcelv/"+river[point]+"/"+pname[point]+"_fft_freq.bin")
    plt.savefig(pm.out_dir()+"/figures/timeseries/"+river[point]+"/"+pname[point]+"_sfcelv.png",dpi=500)
    print "/figures/timeseries/"+river[point]+"/"+pname[point]+"_sfcelv.png"
    #plt.show()
    return 0

#p=Pool(12)
#p.map(make_fig,np.arange(pnum))
#p.terminate()
#--
map(make_fig,np.arange(pnum))
#sendslack(" _finishes_ ")
nc_sfcelv.close()
nc_rmdtrend.close()
nc_rmdseson.close()
nc_standz.close()
