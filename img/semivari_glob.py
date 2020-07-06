#!/opt/local/bin/python
# -*- coding: utf-8 -*-
""" make plots for fitted semivariograms
for each upstream and downstream
Menaka@IIS 2019/10/24"""
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
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#--
def nextXY(x,y,i,j):
  nextxy = pm.CaMa_dir()+"/map/glb_15min/nextxy.bin"
  nextxy = np.fromfile(nextxy,np.int32).reshape(2,720,1440)
  nextx  = nextxy[0]
  nexty  = nextxy[1]
  #--
  lx=[i]
  ly=[j]
  ix=i
  iy=j
  #print i,j
  while (ix!=x or iy!=y):
    iix=ix
    iiy=iy
    ix=nextx[iiy-1,iix-1]
    iy=nexty[iiy-1,iix-1] 
    if (ix==-9 or iy==-9):
      print (ix,iy)
      break  
    lx.append(ix)
    ly.append(iy)
  return lx,ly 
#--
def mk_svfig(sv,lag,pathname,up=0,model="gaussian"):
  # up = 0 is downsream
  # up > 0 is upstream number
  models=['spherical', 'cubic',  'gaussian', 'pentaspherical', 'sineholeeffect','exponential']
  label=['spherical', 'cubic',  'gaussian', 'pentaspherical', 'sineholeeffect','exponential']
  fit=[0,1,1]
  colors=cm.jet(np.linspace(0,1,len(models)))
  #---
  if up == 0:
    iup="downstream"
  else:
    iup="upstream"
  #----
  plt.clf()
  plt.close()
  #fig=plt.figure(figsize=(8.27,11.69))
  fig=plt.figure()#figsize=(11.69,8.27))
  G = gridspec.GridSpec(1,1)
  ax=fig.add_subplot(G[0,0])
  ax.plot(lag,sv,markeredgecolor="black",zorder=101,marker="o",markerfacecolor="none",linestyle="none")#label="downstream",
  N=float(len(lag))
  #--
  #for ii,model in enumerate(models):
  if len(lag)<3:
    a=len(lag)
    ax.plot(lag,sv,color="dodgerblue")
  else: 
#    plsq = LMA.mrqfit(model,lag,sv,fit)
#    ax.plot(lag,LMA.pval(model,lag,plsq[0],fit),color="dodgerblue")#label=label[ii])
    try: 
      plsq = LMA.mrqfit(model,lag,sv,fit)
      ax.plot(lag,LMA.pval(model,lag,plsq[0],fit),color="dodgerblue")#label=label[ii])
      ax.axvline(plsq[0][0],color="dodgerblue",linestyle="--")
      AIC=LMA.AIC_val(model,plsq[0],N,lag,sv,fit)
      aic="%4.2f,%4.2f"%(AIC,plsq[0][0])
      ax.text(0.8,0.02,aic,ha="left",va="center",fontsize=10,transform=ax.transAxes,color="k")
      if lag[-1]<plsq[0][0]:
        a=sum((np.array(lag)<=plsq[0][0])*(1))
      elif plsq[0][0]<0.0:
        a=1
      else:
        a=len(lag)
    except Exception as e:
        print ("exception:",str(e))
        a=1#(np.array(lag)==lag[0])*(1) 
        ax.plot(lag,sv,color="dodgerblue")
  #--title--
  if up <=25:
    stitle= "%s) %s %02d"%(string.ascii_lowercase[up],iup,up)
  else:
    stitle= "%d.%s) %s %02d"%(int(up/26.),string.ascii_lowercase[(up%26)],iup,up)  
  ax.set_title(stitle,fontsize=9,loc="left")
  if int(a)<len(lag):
    ax.set_xlim(xmin=0.0)
  elif len(lag)>2000.0:
    ax.set_xlim(xmin=0.0,xmax=lag[-1]+1.)
  else:
    ax.set_xlim(xmin=0.0,xmax=3000.)
  ax.set_ylim(ymin=0.0,ymax=1.2) 
  #plt.legend(loc="lower right",ncol=1,prop={"size":8})  
  iiup="%02d"%(up)
  figname=pathname+"/"+iup+iiup+".png"
  #plt.show() 
  plt.savefig(figname,dpi=500)   
  return 0#int(a) 
#----------------
def covert_xy_to_latlon(xlist,ylist):
  lat_list=[]  
  lon_list=[]
  for x,y in zip(xlist,ylist):
#    print x, y
    lat_list.append(-(y-1)*0.25 + 90.)
    lon_list.append((x-1)*0.25 - 180.)
  
  return lon_list, lat_list
#-----------------
#def iixy2ixy(lx,ly):
#  low=0
#  high=0
#  for ix in lx:
#    if ix==1:
#      low+=1
#    elif ix==1440:
#      high+=1
#  if low==0 and high==0:
#    return lx,ly
#  elif low<high    
#--
#slink("/hydro/covariance/CaMa_simulation/params.py", "params.py")
#slink("/hydro/covariance/img_code/read_grdc.py","read_grdc.py")
import params as pm
import read_grdc as grdc
#--------------
mk_dir(pm.out_dir()+"/figures")
mk_dir(pm.out_dir()+"/figures/semivar")
#--------------
###nextxy = pm.CaMa_dir()+"/map/glb_15min/nextxy.bin"
###rivwth = pm.CaMa_dir()+"/map/glb_15min/rivwth.bin"
###rivhgt = pm.CaMa_dir()+"/map/glb_15min/rivhgt.bin"
###rivlen = pm.CaMa_dir()+"/map/glb_15min/rivlen.bin"
###elevtn = pm.CaMa_dir()+"/map/glb_15min/elevtn.bin"
###nextxy = np.fromfile(nextxy,np.int32).reshape(2,720,1440)
###rivwth = np.fromfile(rivwth,np.float32).reshape(720,1440)
###rivhgt = np.fromfile(rivhgt,np.float32).reshape(720,1440)
###rivlen = np.fromfile(rivlen,np.float32).reshape(720,1440)
###elevtn = np.fromfile(elevtn,np.float32).reshape(720,1440)
####----
####----
###nextx  = nextxy[0]#ma.masked_where(rivwth<=500.0,nextxy[0]).filled(0)
###nexty  = nextxy[1]#ma.masked_where(rivwth<=500.0,nextxy[1]).filled(0) 
####--------------
models=['spherical', 'cubic',  'gaussian', 'pentaspherical', 'sineholeeffect','exponential']
label=['spherical', 'cubic',  'gaussian', 'pentaspherical', 'sineholeeffect','exponential']
fit=[0,1,1]
colors=cm.jet(np.linspace(0,1,len(models)))
lllat, urlat, lllon, urlon = -90.,90.,-180.,180.
#--
land="#FFFFFF"#"grey"#
water="#C0C0C0"#"royalblue"#
#-- meridians and parallels
meridians = 10.#5.0
parallels = 10.#5.0
#----
pname=[]
xlist=[]
ylist=[]
river=[]
staid=[]
#--
rivernames  = ["AMAZON"]#"LENA"] # ["OB"]#"AMAZONAS"]#["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG"]#,"AMAZONAS","INDUS"] ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS","INDUS"]
#rivernames = grdc.grdc_river_name()
for rivername in rivernames:
#for rivername in ["AMAZONAS"]:#"LENA","NIGER","INDUS","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS"]:
    path = pm.out_dir()+"/figures/semivar/%s"%(rivername)
    print path
    mk_dir(path)
    grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername)
    print rivername, grdc_id,station_loc
    river.append([rivername]*len(station_loc))
    staid.append(grdc_id)
    pname.append(station_loc)
    xlist.append(x_list)
    ylist.append(y_list)


#  oname = "../assim_out/img/sfcelv/%s"%(rivername)
#  mk_dir(oname)
#  station_loc,x_list,y_list = grdc.get_grdc_loc(rivername,"b")
#  for station in station_loc:
#     gid=grdc.get_id(station)
#     if gid== -9999:
#         continue 
#     ix, iy = grdc.get_loc_v394(gid)
#     river.append([rivername])
#     pname.append([station])
#     xlist.append([ix])
#     ylist.append([iy])
#  print rivername, station_loc,x_list,y_list
#  river.append([rivername]*len(station_loc))
#  pname.append(station_loc)
#  xlist.append(x_list)
#  ylist.append(y_list)
#--
#  if rivername=="LENA":
#      river.append([rivername]*3)
#      pname.append(["L1","L2","L3"])
#      xlist.append([1234,1219,1207])
#      ylist.append([  72,  99, 118])
#  if rivername=="NIGER":
#      river.append([rivername]*6)
#      pname.append(["N1","N2","N3","N4","N5","N6"])
#      xlist.append([745,745,733,713,705,701])
#      ylist.append([343,325,311,293,296,304])
#  if rivername=="AMAZONAS":
#      river.append([rivername]*4)
#      pname.append(["B","E","F","G"])
#      xlist.append([516,448,465,421])
#      ylist.append([365,367,417,368])
#  if rivername=="MEKONG":
#      river.append([rivername]*3)
#      pname.append(["Me_D","Me_M","Me_U"])
#      xlist.append([1144,1140,1128])
#      ylist.append([ 320, 294, 283])
#  if rivername=="MISSISSIPPI":
#      river.append([rivername]*4)
#      pname.append(["Mi_D","Mi_M","Mi_U","Mi_U2"]) #最後はミズーリ川
#      xlist.append([362,363,346,309])
#      ylist.append([242,215,138,168])
#  if rivername=="OB":
#      river.append([rivername]*3)
#      pname.append(["O_D","O_M","O_U"])
#      xlist.append([996,997,1049])
#      ylist.append([ 93,122, 160])
#  if rivername=="CONGO":
#      river.append([rivername]*3)
#      pname.append(["C_D","C_M","C_U"])
#      xlist.append([773,814,835])
#      ylist.append([384,354,398])
#  if rivername=="INDUS":
#      river.append([rivername]*3)
#      pname.append(["I_D","I_Sub","I_M"])
#      xlist.append([993,996,1004])
#      ylist.append([252,253, 234])

river=([flatten for inner in river for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])
#***********
print (pm.out_dir())
#-----------
fname=pm.out_dir()+"/semivar/lonlat_list.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
#---
updn={}
#---
for line in lines[1::]:
  line  = filter(None, re.split(" ",line))
  #print line
  lon   = int(line[0]) 
  lat   = int(line[1]) 
  up    = int(line[2])
  dn    = int(line[3])
  #--
  updn[lon,lat]=[up,dn]
#--
pnum=len(pname)
#for point in np.arange(pnum):
def semivari_fig(point):
     lon=xlist[point]+1
     lat=ylist[point]+1
     #--
     up,dn=updn[lon,lat]
     print (pname[point],up,dn)
     #--
     llat=-(lat-1)*0.25 + 90.  
     llon=(lon-1)*0.25 - 180.
     lllat, urlat, lllon, urlon = max(llat-50,-90.),min(llat+50,90.),max(llon-50,-180.),min(llon+50,180.)
     #--
     al=np.zeros([720,1440],np.float32)
     #--
     cs=cm.viridis(np.linspace(0,1,up+2))
     xfig  = int(math.ceil(up/2.0)) + 1
     print (lon,lat)
     #print xfig, cs[0]
     pathname=pm.out_dir()+"/figures/semivar/%s/%s"%(river[point],pname[point])
     #--
     if up + dn == 0:
       #continue
       return 0
     #--
     if dn>0:
       print ("downstream",dn) 
       #--
       lix=[]
       liy=[]
       ldis=[]
       lgamma=[]
       lstd=[] 
       #--
       fname="%s/semivar/%04d%04d/dn%05d.svg"%(pm.out_dir(),lon,lat,0)
       f=open(fname,"r")
       lines=f.readlines()
       f.close()
       #-- 
       for line in lines[1::]:
         line  = filter(None, re.split(" ",line))
         #print line
         ix    = int(line[0])
         iy    = int(line[1]) 
         dis   = float(line[2])
         gamma = float(line[3])
         std   = float(line[4])
         #--
         print (dis, gamma)
         #--
         lix.append(ix)
         liy.append(iy)
         ldis.append(dis)
         lgamma.append(gamma)
         lstd.append(std)
       #---
       #print ldis,lgamma
       #-- 
       mk_dir(pathname)
       mk_svfig(lgamma,ldis,pathname,0)
     #==================================== 
     # upstream 
     if up > 0:
       for iup in np.arange(1,up):
         print ("upstream",iup)
         #--
         lix=[]
         liy=[]
         ldis=[]
         lgamma=[]
         lstd=[] 
         #--
         fname="%s/semivar/%04d%04d/up%05d.svg"%(pm.out_dir(),lon,lat,iup)
         f=open(fname,"r")
         lines=f.readlines()
         f.close()
         #--
         for line in lines[1::]:
           line  = filter(None, re.split(" ",line))
           ix    = int(line[0])
           iy    = int(line[1]) 
           dis   = float(line[2])
           gamma = float(line[3])
           std   = float(line[4])
           #--
           print (dis, gamma)
           #--
           lix.append(ix)
           liy.append(iy)
           ldis.append(dis)
           lgamma.append(gamma)
           lstd.append(std)
         #--
         #print ldis,lgamma
         #--
         mk_svfig(lgamma,ldis,pathname,iup)
     return 0

#print ("parallel")
#p=Pool(5)
#p.map(semivari_fig,np.arange(pnum))
#p.terminate()
#----
map(semivari_fig,np.arange(pnum))

