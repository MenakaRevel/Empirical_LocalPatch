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
import re
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
def grdc_Q(iid,start_dt,last_dt):
  #iid=grdc_id()[name]
  #iid="%d"%(idt)
  # read grdc q
  #grdc ="/cluster/data6/menaka/GRDC_Q/"+iid+".day"
  grdc ="/cluster/data6/menaka/GRDC_2019/"+iid+"_Q_Day.Cmd.txt"
  #--
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  dis = {}
  for line in lines[55::]:
    line     = filter(None, re.split(";",line))
    yyyymmdd = filter(None, re.split("-",line[0]))
    #print yyyymmdd
    yyyy     = int(yyyymmdd[0])
    mm       = int(yyyymmdd[1])
    dd       = int(yyyymmdd[2])
    
    #print start_dt.year,start_dt.month,start_dt.day 
    if start_dt <= datetime.date(yyyy,mm,dd) and last_dt  >= datetime.date(yyyy,mm,dd):
      dis[yyyy,mm,dd]=float(line[2])
      #print float(line[2])
    elif last_dt  < datetime.date(yyyy,mm,dd):
      break

  #---
  start=0
  last=(last_dt-start_dt).days + 1
  Q=[]
  for day in np.arange(start,last):
    target_dt=start_dt+datetime.timedelta(days=day)
    if (target_dt.year,target_dt.month,target_dt.day) in dis.keys():
      Q.append(dis[target_dt.year,target_dt.month,target_dt.day])
    else:
      Q.append(-99.0)
  return np.array(Q)
#--
def grdc_id():
  # read grdc id
  grdc = "/cluster/data6/menaka/GRDC_Q/grdc_list.txt"
  #--
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  #--
  grdc_id={}
  for line in lines[1::]:
    line    = filter(None, re.split(" ",line))
    num     = line[0]
    river   = line[1]
    station = line[3]
    d_info  = line[2]
    print station, num
    grdc_id[station]=num
  return grdc_id
#----
#slink("/hydro/covariance/CaMa_simulation/params.py", "params.py")
#slink("../CaMa_simulation/params.py", "params.py")
#slink("../read_grdc.py","read_grdc.py")
import params as pm
import read_grdc as grdc
#--
mk_dir(pm.out_dir()+"/figures")
mk_dir(pm.out_dir()"/figures/disgraph")
#--read outflow netCDF4--
tag="%04d-%04d"%(pm.starttime()[0],pm.endtime()[0])
fname=pm.out_dir()+"/CaMa_out/"+pm.input_name()+"/outflw"+tag+".nc"
nc=xr.open_dataset(fname)

print pm.patch_start()
syear=pm.starttime()[0]
smonth=pm.starttime()[1]
sdate=pm.starttime()[2]
start_dt=datetime.date(syear,smonth,sdate)
size=60

start=0
#last_dt=datetime.date(int(argvs[1]),int(argvs[2]),int(argvs[3]))
eyear,emonth,edate=pm.patch_end()
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
#rivernames = ["AMAZON"]
rivernames = grdc.grdc_river_name()
for rivername in rivernames:
    path = pm.out_dir()+"/figures/disgraph/%s"%(rivername)
    print path
    mk_dir(path)
    grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername)
    print rivername, grdc_id,station_loc
    river.append([rivername]*len(station_loc))
    staid.append(grdc_id)
    pname.append(station_loc)
    xlist.append(x_list)
    ylist.append(y_list)



###  station_loc,x_list,y_list = grdc.get_grdc_loc(rivername,"b")
###  print rivername, station_loc
###  for station in station_loc:
###    gid=grdc.get_id(station)
###    if gid== -9999:
###        continue
###    #path = "../assim_out/img/disgraph/%s"%(rivername)
###    #print path
###    #mk_dir(path)
###    ix, iy = grdc.get_loc_v394(gid)
###    print station,gid, ix ,iy
###    river.append(rivername)
###    pname.append(station)
###    xlist.append(ix)
###    ylist.append(iy)


#  river.append([rivername]*len(station_loc))
#  pname.append(station_loc)
#  xlist.append(x_list)
#  ylist.append(y_list)
##--
river=([flatten for inner in river for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])
staid=([flatten for inner in staid for flatten in inner])

pnum=len(pname)

#org=[]
#opn=[]
#asm=[]

#swt={}
#for point in np.arange(pnum):
#    swt[point] = []
#
#for day in np.arange(start,last):
#    target_dt=start_dt+datetime.timedelta(days=day)
#    yyyy='%04d' % (target_dt.year)
#    mm='%02d' % (target_dt.month)
#    dd='%02d' % (target_dt.day)
#    #print yyyy,mm,dd
#
#    # make org
#    fname=pm.out_dir()+"/CaMa_out/"+yyyy+"/"+yyyy+mm+dd+"/rivout"+yyyy+".bin"
#    #fname="/media/menaka/HDCL-UT/covariance/CaMa_out/"+yyyy+"/"+yyyy+mm+dd+"/rivout"+yyyy+".bin"
#    orgfile=np.fromfile(fname,np.float32).reshape([720,1440])
#
#    org_frag=[]
#    for point in np.arange(pnum):
#        xpoint=xlist[point]
#        ypoint=ylist[point]
#        org_frag.append(orgfile[ypoint,xpoint])
#        #print orgfile[ypoint,xpoint]
#    org.append(org_frag)
#
#org=np.array(org)
#print np.shape(org),org

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
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    o=np.compress(o>0.0,o)
    s=np.compress(o>0.0,s) 
    return 1 - sum((s-o)**2)/(sum((o-np.mean(o))**2)+1e-20)

#for point in np.arange(pnum):
def make_fig(point):
    plt.close()

    #print org[:,point]
    fig=plt.figure()#figsize=(8.27,11.69))
    G = gridspec.GridSpec(1,1)#4,2)
    ax= fig.add_subplot(G[0,0]) 
    #fig, ax = plt.subplots()
    #org_Q=grdc.grdc_Q(pname[point],start_dt,last_dt)
    #org_Q=np.array(org_Q)

    org_Q=grdc_Q(staid[point],start_dt,last_dt)
    #print org, org_Q
    ed=np.shape(org_Q)[0]
    #print ed , np.shape(org[:,point])
    org=nc.outflw[:,ylist[point],xlist[point]]
    print ed , np.shape(org)
    ax.plot(np.arange(start,last),org,label="CaMa-Flood",color="blue",linewidth=0.7,zorder=102)
    if ed == 0:
        print "no GRDC data"
        #return 0
    else:
        ax.plot(np.arange(start,last),ma.masked_less_equal(org_Q,0.0),label="GRDC",color="black",linewidth=0.7,zorder=101)
        NS1=NS(org,org_Q)
        #NS2=1-((np.sum((org[:ed,point]-org_Q)**2))/(np.sum((org_Q-np.mean(org_Q))**2)))
        #print point,NS1,NS2
        Nash="NS:%4.2f"%(NS1)
        ax.text(0.02,0.95,Nash,ha="left",va="center",transform=ax.transAxes,fontsize=10) 
        plt.legend(loc=1,ncol=1,prop={"size":8})
    # Make the y-axis label, ticks and tick labels match the line color.
    ax.set_ylabel('$discharge$ (m$^3$/s)', color='k')
    ax.set_xlim(xmin=0,xmax=ed+1)
    ax.set_ylim(ymin=0)#,xmax=ed+1)
    ax.tick_params('y', colors='k')
    days=np.arange(start_dt.year,last_dt.year+1,5)
    xxlist=np.linspace(1,N,len(days))
    ax.set_xticks(xxlist)
    ax.set_xticklabels(days,fontsize=8)
    #ax1.set_yticklabels(fontsize=11)
    ax.set_xlabel('$year$', color='k')
    # scentific notaion 
    ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))

    plt.tight_layout(pad=0.2,w_pad=0.05,h_pad=0.05)
    plt.savefig(pm.out_dir()+"/figures/disgraph/"+river[point]+"/"+pname[point]+"_disgraph_GRDC.png",dpi=500)
    #plt.show()
    return 0

#p=Pool(6)
#p.map(make_fig,np.arange(pnum))
#p.terminate()
#--
map(make_fig,np.arange(pnum))
nc.close()
