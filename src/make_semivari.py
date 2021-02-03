#!/opt/local/bin/python
# -*- coding: utf-8 -*-
#==============================================
# Make direcotries for semvari
# Menaka@IIS
# 2020/06/02
import numpy as np
import datetime
import sys
import os
import string
import calendar
import errno
from numpy import ma
import math
import re
#==============================================
def slink(src,dst):
  try:
    os.symlink(src,dst)
  except OSError as exc:  # Python >2.5
    if exc.errno == errno.EEXIST:
      os.remove(dst)
      os.symlink(src,dst)
    else:
      raise
#==============================================
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#==============================================
CaMa_dir=sys.argv[1]
mapname=sys.argv[2]
inname=sys.argv[3]
out_dir=sys.argv[4]
#==============================================
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
f=open(fname,"r")
lines=f.readlines()
f.close()
#==============================================
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
#==============================================
nextxy = CaMa_dir+"/map/"+mapname+"/nextxy.bin"
rivwth = CaMa_dir+"/map/"+mapname+"/rivwth.bin"
rivhgt = CaMa_dir+"/map/"+mapname+"/rivhgt.bin"
rivlen = CaMa_dir+"/map/"+mapname+"/rivlen.bin"
elevtn = CaMa_dir+"/map/"+mapname+"/elevtn.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
#==============================================
nextx  = nextxy[0]
nexty  = nextxy[1]
#==============================================
# make semivar
oname="%s/semivar"%(out_dir)
if not os.path.exists(oname):
  print oname
  mk_dir(oname)
#==============================================
##make sub directory
# /semivar/{mapname}_{inputname}
oname="%s/semivar/%s_%s"%(out_dir,mapname,inname)
if not os.path.exists(oname):
  print oname
  mk_dir(oname)
#==============================================
for i in np.arange(1,nx+1):
  for j in np.arange(1,ny+1):
    if nextx[j-1,i-1]==-9999:
      continue
    oname="%s/semivar/%s_%s/%04d%04d"%(out_dir,mapname,inname,i,j)
    if not os.path.exists(oname):
      print oname
      mk_dir(oname)
#==============================================