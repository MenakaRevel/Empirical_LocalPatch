#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import numpy as np
# import xarray as xr
import re
import os

import params as pm
import my_colorbar as mbar
#=====================
#====================================================================
def vec_par(LEVEL,ax=None):
    ax=ax or plt.gca()
    txt=prename+"_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec "+prename+".txt 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    #print LEVEL, width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    f = open(txt,"r")
    lines = f.readlines()
    f.close()
    #print LEVEL, width, lines, txt
    #---
    for line in lines:
        line    = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        #- higher resolution data
        ixx1 = int((lon1  - west)*60.0)
        iyy1 = int((-lat1 + north)*60.0)

        #----
        ix =catmxy[0,iyy1,ixx1]- 1
        iy =catmxy[1,iyy1,ixx1]- 1

        if ix < 1:
            continue
        
        if rivermap[iy,ix] == 0:
            continue

        if lon1-lon2 > 180.0:
            print (lon1, lon2)
            lon2=180.0
        elif lon2-lon1> 180.0:
            print (lon1,lon2)
            lon2=-180.0
        #--------
        # print data[iy,ix], math.log10(data[iy,ix])
        # colorVal="xkcd:azure"
        #print (lon1,lon2,lat1,lat2,width)
        colorVal=cmap(norm(data[iy,ix]))
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#====================================================================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None,alpha=1):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=110,alpha=alpha,transform=ccrs.PlateCarree())
#====================================================================
# main code
#========================
syear=1979
eyear=2019
inputname="conus_06min_VIC_BC"
outdir="../"
#========================
fname=pm.CaMa_dir()+"/map/"+pm.map_name()+"/params.txt"
with open(fname,"r") as f:
  lines=f.readlines()
#-- map params --
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
west   = float(filter(None, re.split(" ",lines[4]))[0])
east   = float(filter(None, re.split(" ",lines[5]))[0])
south  = float(filter(None, re.split(" ",lines[6]))[0])
north  = float(filter(None, re.split(" ",lines[7]))[0])
#----
nextxy = pm.CaMa_dir()+"/map/"+pm.map_name()+"/nextxy.bin"
rivwth = pm.CaMa_dir()+"/map/"+pm.map_name()+"/rivwth_gwdlr.bin"
rivhgt = pm.CaMa_dir()+"/map/"+pm.map_name()+"/rivhgt.bin"
rivlen = pm.CaMa_dir()+"/map/"+pm.map_name()+"/rivlen.bin"
elevtn = pm.CaMa_dir()+"/map/"+pm.map_name()+"/elevtn.bin"
lonlat = pm.CaMa_dir()+"/map/"+pm.map_name()+"/lonlat.bin"
uparea = pm.CaMa_dir()+"/map/"+pm.map_name()+"/uparea.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
# rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
# rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
# rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
# elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
# uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
rivermap=((nextxy[0]>0))*1.0
#====================================================================
#higher resolution data
catmxy = pm.CaMa_dir()+"/map/"+pm.map_name()+"/1min/1min.catmxy.bin"
catmxy = np.fromfile(catmxy,np.int16).reshape(2,2100,4200)
ix=355
iy=229
loc="%04d%04d"%(ix,iy)
prename="weightage"
# #--read outflow netCDF4--
# tag="%04d-%04d"%(syear,eyear)
# 
# # sfcelv
# fname=outdir+"CaMa_out/"+inputname+"/outflw"+tag+".nc"
# # fname=outdir+"CaMa_out/"+inputname+"/standardized"+tag+".nc"
# nc=xr.open_dataset(fname)

# # projection = ccrs.LambertConformal(central_longitude=-95, central_latitude=40)
# # projection = ccrs.Robinson()

# # f, ax = plt.subplots(subplot_kw=dict(projection=projection))

# # # nc.outflw.mean(dim='time').plot(x="lon",y="lat",ax=ax,transform=ccrs.PlateCarree(), cbar_kwargs=dict(shrink=0.7))
# # nc.outflw.mean(dim='time').plot.imshow(ax=ax,rgb='band',transform=ccrs.PlateCarree(), cbar_kwargs=dict(shrink=0.7))
# data=nc.outflw.mean(dim='time').values
# # data=nc.standardize.mean(dim='time').values
fname="../weightage/"+inputname+"_1000KM_dam/"+loc+".bin"
data=np.fromfile(fname,np.int32).reshape(ny,nx)
rivermap=rivermap*(data>1e-20)
data=data*rivermap
# data=np.ma.masked_less(data,0.00)

land="white"
water="white"

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

# cmap=cm.get_cmap('RdYlBu')
cmap=mbar.colormap("H01")
# cmap=cm.get_cmap("Blues")
# vmin=1.0
# vmax=2.0e5
# norm=LogNorm(vmin=vmin,vmax=vmax)

vmin=0.0
vmax=1.0
norm=Normalize(vmin=vmin,vmax=vmax)

# river width
sup=2
w=0.01
alpha=1
width=0.5
#------
print (west,east,south,north)
resol=1
fig=plt.figure()
G = gridspec.GridSpec(1,1)
# ax=fig.add_subplot(G[0,0])
# ax = fig.add_subplot(G[0,0],projection=ccrs.Robinson()) 
ax = fig.add_subplot(G[0,0],projection=ccrs.LambertConformal())
#-----------------------------  
ax.set_extent([west,east,south,north],crs=ccrs.PlateCarree())
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', edgecolor="k", facecolor=land, linewidth=0.1),zorder=100)
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+pm.CaMa_dir()+" "+pm.map_name()+" > "+prename+".txt")  
map(vec_par,np.arange(1,10+1,1))
# map(vec_par,np.arange(6,10+1,1))
# map(vec_par,np.arange(1,5+1,1))
# ix , iy
lon=lonlat[0,iy-1,ix-1]
lat=lonlat[1,iy-1,ix-1]
print (lon, lat)
ax.scatter(lon,lat,s=20,marker="o",color="green",transform=ccrs.PlateCarree(),zorder=105)
# ax.outline_patch.set_linewidth(0.0)

# ax.coastlines()
im=plt.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm,zorder=101)
im.set_visible(False)
# cax=fig.add_axes([0.4,0.25,0.35,0.01])
l,b,w,h=ax.get_position().bounds
cax=fig.add_axes([l,b-0.01*h,0.8*w,0.01])
cbar=plt.colorbar(im,extend='max',orientation="horizontal",cax=cax)
cbar.ax.tick_params(labelsize=6)
cbar.set_label("Spatial Depndancey Weight",fontsize=8)

plt.show()
# plt.savefig(outdir+"CaMa_out/"+inputname+"/weightage_"+loc+".jpg")
os.system("rm -r "+prename+"*.txt")