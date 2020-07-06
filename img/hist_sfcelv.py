#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import sys
import os
import errno
import math
from scipy import stats
from numpy import ma 
import matplotlib.gridspec as gridspec
import string
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
#---
#--
#slink("/hydro/covariance/CaMa_simulation/params.py", "params.py")
#slink("/hydro/SWOTDA/img_code/read_grdc.py","read_grdc.py")
import params as pm
import read_grdc as grdc
#--
#--
mk_dir(pm.out_dir()+"/figures")
mk_dir(pm.out_dir()+"/figures/hist_sfcelv")
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

#colors = pc.color_assimd()
pname=[]
xlist=[]
ylist=[]
river=[]
staid=[]
#--
#rivernames  = ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS","INDUS"] 
rivernames=["AMAZON"]#
#rivernames = ["MEKONG"]
#rivernames = grdc.grdc_river_name()
for rivername in rivernames:
    path = pm.out_dir()+"/figures/hist_sfcelv/%s"%(rivername)
    print path
    mk_dir(path)
    grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername)
    print rivername, grdc_id,station_loc
    river.append([rivername]*len(station_loc))
    staid.append(grdc_id)
    pname.append(station_loc)
    xlist.append(x_list)
    ylist.append(y_list)

  ##path = "../img/hist_sfcelv/%s"%(rivername)
  ##print path
  ##mk_dir(path) 
  ##station_loc,x_list,y_list = grdc.get_grdc_loc(rivername,"b")
  ##print rivername, station_loc
  ##river.append([rivername]*len(station_loc))
  ##pname.append(station_loc)
  ##xlist.append(x_list)
  ##ylist.append(y_list)
#--
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

###org=[]
###tdl=[]
###fft=[]
###sdz=[]
###for day in np.arange(start,last):
###    target_dt=start_dt+datetime.timedelta(days=day)
###    yyyy='%04d' % (target_dt.year)
###    mm='%02d' % (target_dt.month)
###    dd='%02d' % (target_dt.day)
###    print yyyy,mm,dd
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
###print np.shape(org)

def NS(s,o):
    """
    Nash Sutcliffe efficiency coefficient
    input:
        s: simulated
        o: observed
    output:
        ns: Nash Sutcliffe efficient coefficient
    """
    #s,o = filter_nan(s,o)
    o=ma.maksed_where(o<=0.0,o).filled(0.0)
    s=ma.maksed_where(o<=0.0,s).filled(0.0) 
    return 1 - sum((s-o)**2)/sum((o-np.mean(o))**2)
#------------------
def KolmogorovSmirnovVal(N,a):
    # N number of data
    # a alpha [0.1,0.05,0.02,0.01]
    #--
    kst={0.1:1.22,0.05:1.36,0.02:1.51,0.01:1.63}
    
    CV=kst[a]/math.sqrt(N)
    return CV
#------------------    

#for point in np.arange(10):#pnum):

def make_fig(point):
    plt.close()

    #print org[:,point]
    # standardized data
    sdz=nc_standz.standardize[:,ylist[point],xlist[point]]
    # water surface elevation
    org=nc_sfcelv.sfcelv[:,ylist[point],xlist[point]]
    # seasonality removed
    fft=nc_rmdseson.rmdsesn[:,ylist[point],xlist[point]]
    # trend line
    tdl=nc_rmdtrend.rmdtrnd[:,ylist[point],xlist[point]]
    #fig=plt.figure()
    fig=plt.figure(figsize=(8.27,11.69))
    G = gridspec.GridSpec(3,1)
    ax1 = fig.add_subplot(G[0,0])
    mu = np.mean(org.values)
    sigma=np.std(org.values)
    bin1=100#int(abs(np.amax(org[:,point])-np.amin(org[:,point])))
    print mu, sigma, bin1#, type(mu)
    #print np.amax(org[:,point]),np.amin(org[:,point])
    n, bins, patches = ax1.hist(org, bin1, normed=1, facecolor='green', alpha=0.75)#-np.mean(org[:,point]
    # add a 'best fit' line
    y = mlab.normpdf( bins, mu, sigma)
    l = ax1.plot(bins, y, 'r--', linewidth=1)
    ms="mean:%6.2f\nstd:%6.2f"%(mu,sigma)
    ax1.text(0.02,0.9,ms,ha="left",va="center",fontsize=8,transform=ax1.transAxes,color="black")
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('$frequency$', color='k')
    ax1.set_xlabel('WSE ($m$)', color='k')
    #ax1.set_xlim(xmin=0,xmax=la+1)
    #ax1.set_ylim(ymin=0)#,xmax=ed+1)
    ax1.tick_params('y', colors='k')
    #------------------------
    # seasonality removed data 
    ax2 = fig.add_subplot(G[1,0])
#    print np.log(org[:,point])
    mu = np.mean(fft.values)
    sigma=np.std(fft.values)
    bin1=100#int(abs(np.amax(org[:,point])-np.amin(org[:,point])))
    print mu, sigma, bin1
    #print np.amax(org[:,point]),np.amin(org[:,point])
    n, bins, patches = ax2.hist(fft, bin1, normed=1, facecolor='green', alpha=0.75)#-np.mean(org[:,point]
    # add a 'best fit' line
    y = mlab.normpdf( bins, mu, sigma)
    l = ax2.plot(bins, y, 'r--', linewidth=1)
    ms="mean:%6.2f\nstd:%6.2f"%(mu,sigma)
    ax2.text(0.02,0.9,ms,ha="left",va="center",fontsize=8,transform=ax2.transAxes,color="black")
    # Make the y-axis label, ticks and tick labels match the line color.
    ax2.set_ylabel('$frequency$', color='k')
    ax2.set_xlabel('WSE - TD -SEASON', color='k')
    #ax2.set_xlim(xmin=0,xmax=la+1)
    #ax2.set_ylim(ymin=0)#,xmax=ed+1)
    ax2.tick_params('y', colors='k')
    #------------------------
    # standadized data 
    ax3 = fig.add_subplot(G[2,0])
#    print np.log(org[:,point])
    mu = np.mean(sdz.values)
    sigma=np.std(sdz.values)
    #ks=stats.kstest(sdz[:,point],"norm",stats.norm.fit(sdz[:,point]))
    #ks=stats.kstest(sdz[:,point],"norm",args=(mu,sigma))
    #ks=stats.normaltest(sdz[:,point].T)
    ks=stats.shapiro(sdz.values)
    bin1=100#int(abs(np.amax(sdz[:,point])-np.amin(sdz[:,point])))
    print mu, sigma, bin1
    #print np.amax(sdz[:,point]),np.amin(sdz[:,point])
    n, bins, patches = ax3.hist(sdz, bin1, normed=1, facecolor='green', alpha=0.75)#-np.mean(sdz[:,point]
    # add a 'best fit' line
    y = mlab.normpdf( bins, mu, sigma)
    l = ax3.plot(bins, y, 'r--', linewidth=1)
    ms="mean:%6.2f\nstd:%6.2f\nks:%6.4f\np-val:%0.2E"%(mu,sigma,ks[0],ks[1])
    print ms
    ax3.text(0.02,0.7,ms,ha="left",va="center",fontsize=8,transform=ax3.transAxes,color="black")
    # Make the y-axis label, ticks and tick labels match the line color.
    ax3.set_ylabel('$frequency$', color='k')
    ax3.set_xlabel('$standardized$ $data$', color='k')
    #ax3.set_xlim(xmin=0,xmax=la+1)
    #ax3.set_ylim(ymin=0)#,xmax=ed+1)
    ax3.tick_params('y', colors='k')

    
    #ax.text(0.02,0.9,Nash,ha="left",va="center",transform=ax.transAxes,fontsize=8) 
    #plt.legend(loc=1)
    plt.savefig(pm.out_dir()+"/figures/hist_sfcelv/"+river[point]+"/"+pname[point]+"_hist_sfcelv.png",dpi=500)
    #plt.show()
    return 0

#p=Pool(10)
#p.map(make_fig,np.arange(pnum))
#p.terminate()
#---
map(make_fig,np.arange(pnum))

for a in [0.1,0.05,0.02,0.01]:
  print a,KolmogorovSmirnovVal(N,a)

def stanz_fig(point):
    plt.close()

    #print org[:,point]
    # standardized data
    sdz=nc_standz.standardize[:,ylist[point],xlist[point]]
    #fig=plt.figure()
    fig=plt.figure(figsize=(8.27,4))#11.69/3.0))
    G = gridspec.GridSpec(1,1)
    # standadized data 
    ax3 = fig.add_subplot(G[0,0])
#    print np.log(org[:,point])
    mu = np.mean(sdz.values)
    sigma=np.std(sdz.values)
    #ks=stats.kstest(sdz[:,point],"norm",stats.norm.fit(sdz[:,point]))
    #ks=stats.kstest(sdz[:,point],"norm",args=(mu,sigma))
    #ks=stats.normaltest(sdz[:,point].T)
    ks=stats.shapiro(sdz)
    bin1=100#int(abs(np.amax(sdz[:,point])-np.amin(sdz[:,point])))
    print mu, sigma, bin1
    #print np.amax(sdz[:,point]),np.amin(sdz[:,point])
    n, bins, patches = ax3.hist(sdz, bin1, normed=1, facecolor='green', alpha=0.75)#-np.mean(sdz[:,point]
    # add a 'best fit' line
    y = mlab.normpdf( bins, mu, sigma)
    l = ax3.plot(bins, y, 'r--', linewidth=1)
    #ms="mean:%6.2f\nstd:%6.2f\nks:%6.4f\np-val:%0.2E"%(mu,sigma,ks[0],ks[1])
    ms="mean:%6.2f\nstd:%6.2f"%(mu,sigma)
    print ms
    ax3.text(0.02,0.8,ms,ha="left",va="center",fontsize=10,transform=ax3.transAxes,color="black")
    # Make the y-axis label, ticks and tick labels match the line color.
    ax3.set_ylabel('$frequency$', color='k')
    ax3.set_xlabel('$standardized$ $data$', color='k')
    #ax3.set_xlim(xmin=0,xmax=la+1)
    #ax3.set_ylim(ymin=0)#,xmax=ed+1)
    ax3.tick_params('y', colors='k')
    plt.savefig(pm.out_dir()+"/figures/hist_sfcelv/"+river[point]+"/"+pname[point]+"_hist_standz.png",pad_inches=0.02,dpi=800)
    return 0

#p=Pool(10)
#p.map(stanz_fig,np.arange(pnum))
#p.terminate()
#---
map(stanz_fig,np.arange(pnum))
#close the netCDF4 files
nc_rmdtrend.close()
nc_rmdseson.close()
nc_standz.close()
nc_sfcelv.close()
