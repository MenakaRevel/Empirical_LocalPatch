#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib import colors
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
from mpl_toolkits.basemap import Basemap
#from slacker import Slacker
from multiprocessing import Pool
from multiprocessing import Process
#---
os.system("cp ../params.py params.py")
import LMA_semivari as LMA
import params as pm
import read_grdc as grdc
import slice_region as sl
import river_boundry as rbd
#-----------------------------
def latlon_river(rivername,ix,iy,mapname="glb_06min",nYY=1800,nXX=3600):
    global lllat, urlat, lllon, urlon
    lonlat = pm.CaMa_dir()+"/map/glb_"+mapname[-5::]+"/lonlat.bin"
    lonlat = np.fromfile(lonlat,np.float32).reshape(2,nYY,nXX)
    llon=lonlat[0,iy-1,ix-1]
    llat=lonlat[1,iy-1,ix-1]
    adj=20.0
    lllat, urlat, lllon, urlon = max(llat-adj,-90.),min(llat+adj,90.),max(llon-adj,-180.),min(llon+adj,180.)
    if rivername=="LENA":
        lllat = 50.
        urlat = 80.
        lllon = 100.
        urlon = 145.
    if rivername=="NIGER":
        lllat = 0.
        urlat = 25.
        lllon = -10.
        urlon = 15.
    if rivername=="AMAZONAS" or rivername=="AMAZON":
        lllat = -20.
        urlat = 10.
        lllon = -80.
        urlon = -45.
    if rivername=="MEKONG":
        lllat = 10.
        urlat = 35.
        lllon = 90.
        urlon = 120.
    if rivername=="MISSISSIPPI":
        lllat = 20.
        urlat = 50.
        lllon = -115.
        urlon = -75.
    if rivername=="OB":
        lllat = 40.
        urlat = 70.
        lllon = 55.
        urlon = 95.
    if rivername=="CONGO":
        lllat = -15.
        urlat = 10.
        lllon = 10.
        urlon = 35.
    if rivername=="INDUS":
        lllat = 20.
        urlat = 40.
        lllon = 60.
        urlon = 80.
    if rivername=="VOLGA":
        lllat = 40.
        urlat = 65.
        lllon = 30.
        urlon = 70.
    if rivername=="NILE":
        lllat = -5.
        urlat = 30.
        lllon = 20.
        urlon = 40.
    if rivername=="YUKON":
        lllat = 55.
        urlat = 75.
        lllon = -165.
        urlon = -130.
    #if rivername not in ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS","INDUS"]:
    #    adj=20.
    #    lllat, urlat, lllon, urlon = max(llat-adj,-90.),min(llat+adj,90.),max(llon-adj,-180.),min(llon+adj,180.)

    return lllat, urlat, lllon, urlon
#----
def riveridname(rivername):
    river=rivername[0]+rivername[1::].lower()
    if rivername=="LENA":
        river="Lena"
    if rivername=="NIGER":
        river="Niger"
    if rivername=="AMAZONAS":
        river="Amazon"
    if rivername=="MEKONG":
        river="Mekong"  
    if rivername=="MISSISSIPPI":
        river="Mississippi"
    if rivername=="OB":
        river="Ob"
    if rivername=="CONGO":
        river="Congo"
    if rivername=="INDUS":
        river="Indus"
    if rivername=="ST._LAWRENCE":
        river="St Lawrence"
    if rivername=="BRAHMAPUTRA":
        river="Ganges-Brahamaputra"

    return river
#----
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#--
fname=pm.CaMa_dir()+"/map/"+pm.map_name()+"/params.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
#-- map params --
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
west   = float(filter(None, re.split(" ",lines[4]))[0])
east   = float(filter(None, re.split(" ",lines[5]))[0])
south  = float(filter(None, re.split(" ",lines[6]))[0])
north  = float(filter(None, re.split(" ",lines[7]))[0])
# regional map dimesion realtive to global map
offsetX= int((west+180.0)/gsize)
offsetY= int((90.0-north)/gsize)
# global dimension
glbname="glb_"+pm.map_name()[-5::]
nXX    = int((180.0+180.0)/gsize)
nYY    = int((90.0+90.0)/gsize)
print glbname
#--------------
nextxy = pm.CaMa_dir()+"/map/"+pm.map_name()+"/nextxy.bin"
rivwth = pm.CaMa_dir()+"/map/"+pm.map_name()+"/rivwth.bin"
rivhgt = pm.CaMa_dir()+"/map/"+pm.map_name()+"/rivhgt.bin"
rivlen = pm.CaMa_dir()+"/map/"+pm.map_name()+"/rivlen.bin"
elevtn = pm.CaMa_dir()+"/map/"+pm.map_name()+"/elevtn.bin"
uparea = pm.CaMa_dir()+"/map/"+pm.map_name()+"/uparea.bin"
lonlat = pm.CaMa_dir()+"/map/"+pm.map_name()+"/lonlat.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
#----
rivnum = pm.out_dir()+"/dat/rivnum_"+pm.map_name()[-5::]+".bin"
rivnum = np.fromfile(rivnum,np.int32).reshape(nYY,nXX)
#----
nextx  = nextxy[0]#ma.masked_where(rivwth<=500.0,nextxy[0]).filled(0)
nexty  = nextxy[1]#ma.masked_where(rivwth<=500.0,nextxy[1]).filled(0) 
#--------------
#fname = "./congo_nextxy.bin"
#c_nextx= np.fromfile(fname,np.int32).reshape([2,nyYYnXX])[0]
#fname = pm.CaMa_dir()+"/map/glb_15min/outclm.bin"
#trueforo = np.fromfile(fname,np.float32).reshape([2,nYY,nXX])
#dis=(trueforo[0]>500.)*1
#dis=dis*((nextx>0)*1)
#--major rivers and Ids
rivid={}
fname=pm.out_dir()+"/dat/river30_id.txt"
f = open(fname,"r")
lines = f.readlines()
f.close()
#---
for line in lines:
  line    = filter(None, re.split(",",line))
  riverid = int(line[0])
  river   = filter(None, re.split("\n",line[1]))[0].strip()
  #print river 
  rivid[river]=riverid
#----
pname=[]
xlist=[]
ylist=[]
river=[]
staid=[]
#--
#rivernames  = ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS"]#,"INDUS"]# ["LENA","NIGER","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS","INDUS"] ["AMAZONAS"]#"CONGO"]#
#rivernames = grdc.grdc_river_name()
rivernames = ["AMAZON"]
for rivername in rivernames:
#for rivername in ["AMAZONAS"]:#"LENA","NIGER","INDUS","CONGO","OB","MISSISSIPPI","MEKONG","AMAZONAS"]:
#  oname = "../assim_out/img/sfcelv/%s"%(rivername)
#  mk_dir(oname)
    grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername)
    print rivername, grdc_id,station_loc
    river.append([rivername]*len(station_loc))
    staid.append(grdc_id)
    pname.append(station_loc)
    xlist.append(x_list)
    ylist.append(y_list)
###
###    station_loc,x_list,y_list = grdc.get_grdc_loc(rivername,"b")
###    for station in station_loc:
###        gid=grdc.get_id(station)
###        if gid== -9999:
###            continue 
###        ix, iy = grdc.get_loc_v394(gid)
###        river.append([rivername])
###        pname.append([station])
###        xlist.append([ix-offsetX])
###        ylist.append([iy-offsetY])
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
##  if rivername=="CONGO":
##      river.append([rivername]*3)
##      pname.append(["C_D","C_M","C_U"])
##      xlist.append([773,814,835])
##      ylist.append([384,354,398])
##  if rivername=="INDUS":
##      river.append([rivername]*3)
##      pname.append(["I_D","I_Sub","I_M"])
##      xlist.append([993,996,1004])
##      ylist.append([252,253, 234])
#
#  if rivername=="CONGO":
#      river.append([rivername]*6)
#      pname.append(["C1","C2","C3","C4","C5","C6"])
#      xlist.append([834,813,794,795,784,772])
#      ylist.append([397,353,342,373,372,383])


river=([flatten for inner in river for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])
#***********
#print "out"#pm.out_dir()
#---
cmap=colors.ListedColormap(["#C0C0C0","dodgerblue"])
bounds=[-1,0.2,1]
norm=colors.BoundaryNorm(bounds,cmap.N)
cmap1=colors.ListedColormap(["white","#C0C0C0"])
bounds1=[1.,0.0,-1.]
norm1=colors.BoundaryNorm(bounds1,cmap1.N)
#--
land="#FFFFFF"#"grey"#
water="#C0C0C0"#"royalblue"#
#-- meridians and parallels
meridians = 5.0
parallels = 5.0
#--
w=0.15*2
alpha=1
alpha1=1
width=0.5
#---
#---
local_patch="local_patch"#"_%3.2f"%(pm.threshold())
local_patch1="local_patch_one_%3.2f"%(pm.threshold())
#pathname="../img/local_patch_one_0.90"
pathname="../figures/%s"%(local_patch1)
mk_dir(pathname)
pnum=len(pname)
#--
#local_patch="local_patch_0.90"
for point in np.arange(pnum):
#def mk_fig(point):
  ix=xlist[point]+1
  iy=ylist[point]+1
  #--
  if not riveridname(river[point]) in rivid.keys():
    continue
  #----
  fname=pm.out_dir()+"/weightage/%04d%04d.bin"%(ix,iy)
  wgt=np.fromfile(fname,np.float32).reshape([ny,nx])
  #--
  rivername=river[point]
  c_nextx=(rivnum==rivid[riveridname(rivername)])*1.0
  #--
#  x=int(point/2)
#  y=point%2
#  print point,x,y 
  plt.close() 
#--figure in A4 size
  #fig=plt.figure(figsize=(8.27,11.69))
  fig=plt.figure()#figsize=(11.69,8.27))
  G = gridspec.GridSpec(1,1)
  #  mtitle="Assimilation Index"
  #  fig.suptitle(mtitle,fontsize=12,fontweight="bold")
  latlon_river(rivername,ix,iy)
  #print rivername
  #lllat,urlat,lllon,urlon=rbd.latlon_river(rivername)
  ax1 = fig.add_subplot(G[0,0])
  #data=sl.slice(wgt,lllat,urlat,lllon,urlon,gsize)
  #data_Q=sl.slice(dis,lllat,urlat,lllon,urlon,gsize)
  #data=ma.masked_greater_equal(data,0.8).filled(1.0)
  M  = Basemap(resolution="c", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon,ax=ax1) 
#  M.drawcoastlines(color="k",linewidth=0.5, zorder=100)
  M.fillcontinents(color=land,lake_color=water,zorder=99)
  M.drawmapboundary(fill_color=water,zorder=98)
  #--
  M.drawmeridians(np.arange(lllon,urlon+1, meridians), labels=[0, 0, 0, 1], fontsize=14, rotation=0,linewidth=0.5,zorder=102)
  M.drawparallels(np.arange(lllat,urlat+1, parallels), labels=[1, 0, 0, 0], fontsize=14,linewidth=0.5,zorder=102)
#  M.imshow(ma.masked_less_equal(data_Q,0.0),interpolation="nearest",cmap=cmap1,origin="upper",zorder=100,norm=norm1)
#  im=M.imshow(ma.masked_less_equal(data,0.6),interpolation="nearest",cmap=cmap,origin="upper",zorder=101,norm=norm)
  #-----------
  #--
  box="%f %f %f %f"%(lllon,urlon,urlat,lllat)
  #  os.system("./bin/txt_vector "+str(lllon)+str(urlon)+str(urlat)+str(lllat)+" > tmp.txt")
  os.system("./bin/txt_vector "+box+" "+pm.CaMa_dir()+" "+glbname+" > tmp.txt")
  for LEVEL in range(1,10+1):
    os.system("./bin/print_rivvec tmp.txt 1 "+str(LEVEL)+" > tmp2.txt")
    width=float(LEVEL)*w
    #print width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    f = open("tmp2.txt","r")
    lines = f.readlines()
    f.close()
    #----------------
    for line in lines:
        line = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])
        #----
        iix  = int((lon1 + 180.)*(1.0/gsize))
        iiy  = int((-lat1 + 90.)*(1.0/gsize))
        #----
        if c_nextx[iiy,iix] <= 0:
          continue
        #print lon1,lat1,width
        x1,y1=M(lon1,lat1)
        x2,y2=M(lon2,lat2)
        M.plot([x1,x2],[y1,y2],color="#C0C0C0",linewidth=width,zorder=101,alpha=alpha)
  #--
  fname=pm.out_dir()+"/"+local_patch+"/patch%04d%04d.txt"%(ix,iy)
  f=open(fname,"r")
  lines = f.readlines()
  f.close()
  #--
  f=open("tmp.txt","w")
  for line in lines:#[:1]:
    line = filter(None, re.split(" ",line))
    #print line
    iix = int(line[0])
    iiy = int(line[1])
    jjx = nextx[iiy-1,iix-1]
    jjy = nexty[iiy-1,iix-1]
    #-----
    lon1 = lonlat[0,iiy-1,iix-1]
    lat1 = lonlat[1,iiy-1,iix-1]
    lon2 = lonlat[0,jjy-1,jjx-1]
    lat2 = lonlat[1,jjy-1,jjx-1]
    #--
    line1="%12.5f %12.5f %12.5f %12.5f %12.1f\n"%(lon1, lat1, lon2, lat2, uparea[iiy-1,iix-1]/1000.**2)
    f.write(line1)
  f.close()
  for LEVEL in range(1,10+1):
    os.system("./bin/print_rivvec tmp.txt 1 "+str(LEVEL)+" > tmp2.txt")
    width=float(LEVEL)*w
    #print width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    f = open("tmp2.txt","r")
    lines = f.readlines()
    f.close() 
    #--- 
    for line in lines:
        line = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])
        #----
        iix = int((lon1 + 180.)*(1.0/gsize))
        iiy = int((-lat1 + 90.)*(1.0/gsize))
        #----
        if c_nextx[iiy,iix] <= 0:
          continue
        #print lon1,lat1,lon2,lat2,width
        x1,y1=M(lon1,lat1)
        x2,y2=M(lon2,lat2)
        M.plot([x1,x2],[y1,y2],color="dodgerblue",linewidth=width,zorder=102,alpha=alpha)
  # target pixel 
  lon = west  + (ix-1)*gsize
  lat = north - (iy-1)*gsize
  plt.scatter(lon,lat,s=100,marker="o",color="red",zorder=105)
  #plt.annotate(annotate_string,xy=(lon,lat),xycoords="data",horizontalalignment="left",verticalalignment="top",fontsize=12,zorder=106)
  #--title--
#  stitle= "%s)"%(string.ascii_lowercase[point])
#  ax1.set_title(stitle,fontsize=9,loc="left")
  #--
  pathname="../figures/"+local_patch1+"/"+rivername
  #pathname=pathname+"/"+rivername
  mk_dir(pathname)
  figname=pathname+"/"+pname[point]+".png"
  print (figname)
  plt.savefig(figname)
  #plt.show()
# remove tmp
os.system("rm -r tmp*")
