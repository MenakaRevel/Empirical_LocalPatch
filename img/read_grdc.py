# usr/lib/python 
""" read GRDC gauaging locations"""
import numpy as np
import re
import shutil
import os
#--
#os.system("ln -sf ../params.py params.py")
#shutil.copy("../params.py","params.py")
import params as pm
#--
def get_grdc_loc(name,info = "a"):
  #--ask the river name and a or b
  # a - most downsream loc
  # b - all locations  
  grdc = pm.CaMa_dir() + "/map/glb_15min/grdc_loc.txt"
  
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  station_loc = []
  x_list      = []
  y_list      = [] 
  #---
  for line in lines:
    line    = filter(None, re.split(" ",line))
    grdc_id = line[0]
    river   = line[1]
    station = line[2]
    d_info  = line[3]
    #--
    lon     = float(line[4])
    lat     = float(line[5])
    #--
    ix      = int(line[6])-1
    iy      = int(line[7])-1
    
    if name == river:
         if d_info == "a" and info == "a": #info == d_info and info == "a":
             #print d_info
             station_loc.append(station)
             x_list.append(ix)
             y_list.append(iy)
         elif info == "b"  and info in ["a","b"]:
             #print d_info
             station_loc.append(station)
             x_list.append(ix)
             y_list.append(iy)
  
  return station_loc,x_list,y_list
#------------
def get_grdc_loc_dic(name,info = "a"):
  #--aske the river name and a or b
  # a - most downsream loc
  # b - all locations  
  grdc = pm.CaMa_dir() + "/map/glb_15min/grdc_loc.txt"
  
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  station_loc = {}
#  print station_loc 
  for line in lines:
    line    = filter(None, re.split(" ",line))
    grdc_id = line[0]
    river   = line[1]
    station = line[2]
    d_info  = line[3]
    station_loc[station] = ()
    #--
    lon     = float(line[4])
    lat     = float(line[5])
    #--
    ix      = int(line[6])-1
    iy      = int(line[7])-1
    if info == d_info and info == "a":
      station_loc[station] = (ix,iy)
    else:
      station_loc[station] = (ix,iy)
  
#  print station_loc   
  return station_loc   

#----
def grdc_river_name():
  grdc = pm.CaMa_dir() + "/map/glb_15min/grdc_loc.txt"
  
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  rivername =  []
  
  for line in lines:
    line    = filter(None, re.split(" ",line))
    grdc_id = line[0]
    river   = line[1]
    d_info  = line[3]
    if d_info == "a":
      rivername.append(river)

  return rivername
#--
def get_grdc_loc_latlon(name,info = "a"):
  #--aske the river name and a or b
  # a - most downsream loc
  # b - all locations  
  grdc = pm.CaMa_dir() + "/map/glb_15min/grdc_loc.txt"
  
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  station_loc = []
  x_list      = []
  y_list      = [] 
  #---
  for line in lines:
    line    = filter(None, re.split(" ",line))
    grdc_id = line[0]
    river   = line[1]
    station = line[2]
    d_info  = line[3]
    #--
    lon     = float(line[4])
    lat     = float(line[5])
    #--
    ix      = int(line[6])-1
    iy      = int(line[7])-1
    if info == d_info and info == "a":
      station_loc.append(station)
      x_list.append(lon)
      y_list.append(lat)
    else:
      station_loc.append(station)
      x_list.append(lon)
      y_list.append(lat)
  
  return station_loc,x_list,y_list
#---
def get_id(name):
  #--get GRDC id
  grdc = pm.CaMa_dir() + "/map/glb_15min/grdc_loc.txt"
  
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  
  gid=-9999
  #---
  for line in lines:
    line    = filter(None, re.split(" ",line))
    grdc_id = int(line[0])
    river   = line[1]
    station = line[2]
    d_info  = line[3]
    #--
      
    if name == station:
      gid=grdc_id

  return gid
#--
def get_loc_v394(gid):
  #--ask the river name and a or b
  # a - most downsream loc
  # b - all locations  
  grdc = pm.CaMa_dir() + "/map/glb_15min/GRDC_alloc.txt"
  
  f = open(grdc,"r")
  lines = f.readlines()
  f.close()
  ix=0
  iy=0
  gid=int(gid)
  #---
  for line in lines[1::]:
    line    = filter(None, re.split(" ",line))
    grdc_id = int(line[0])
    u_info  = line[7]
    #--
    if gid==grdc_id:
      ix      = int(line[8])-1
      iy      = int(line[9])-1
  
  return ix,iy
#---
def grdc_Q(name,start_dt,last_dt):
  iid=grdc_id()[name]
  #iid="%d"%(idt)
  # read grdc q
  grdc ="/cluster/data6/menaka/covariance/GRDC_Q/"+iid+".day"
  
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
  grdc = "/cluster/data6/menaka/covariance/GRDC_Q/grdc_list.txt"
  
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
    #print station, num
    grdc_id[station]=num
  return grdc_id
