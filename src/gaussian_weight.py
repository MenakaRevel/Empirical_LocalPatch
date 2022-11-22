#!/opt/local/bin/python
# -*- coding: utf-8 -*-
""" fit the semivariogams using gaussian model
  use LMA_semivari to fit the model
  Gaussian weight => w= exp(-r^2 / 2 * sigma^2) see Miyoshi etal 2007a,b
  ******Calculate the gaussian weights only*********
  Menaka@IIS 2021/08/24"""
import numpy as np
import datetime
import sys
import os
import string
import calendar
import errno
import re
import math
from numpy import ma
from multiprocessing import Pool
from multiprocessing import Process
#==============================================
import LMA_semivari as LMA
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
def gaussian(i, p):
    p=abs(p)
    #yfit=[]
    #for i in x:
    i=abs(i)
    if i>=0:
      yfit=math.exp(-3.0*(i/p)**2.0)
    else:
      yfit=0.
    return yfit
#==============================================
def Gauss_w(x,s=500.):
  yfit=[]
  for i in x:
    i=abs(i)
    yfit.append((math.exp((-i**2.0)/(2*s**2.0))))
  return np.array(yfit)
#==============================================
def mk_svg1(sv,lag,model="gaussian"):
  # up = 0 is downsream
  # up > 0 is upstream number
  fit=[0,1,1]
  #----
  N=float(len(lag))
  #--
  if len(lag)<3:
    a=lag[-1]
  else:
    try: 
      plsq = LMA.mrqfit(model,lag,sv,fit)
      if lag[-1]>=plsq[0][0]:
        a=plsq[0][0]
      elif plsq[0][0]<0.:
        a=lag[0]
      else:
        a=lag[-1]
    except Exception as e:
        print "except:",str(e)
        a=lag[0] 
  return a
#==============================================
def mk_svg(sv,lag,model="gaussian"):
  # up = 0 is downsream
  # up > 0 is upstream number
  fit=[0,1,1]
  #----
  N=float(len(lag))
  #--
  if len(lag)<3:
    c,a=sv[-1],lag[-1]
  else:
    try: 
      plsq = LMA.mrqfit(model,lag,sv,fit)
      if lag[-1]>=plsq[0][0] and plsq[0][0]>0.0 and plsq[0][1]>0.0:
        c,a=plsq[0][1],plsq[0][0]
      elif plsq[0][0]<=0.0 or plsq[0][1]<=0.0:
        c,a=sv[-1],lag[-1]
      else:
        c,a=sv[-1],lag[-1]
    except Exception as e:
        print "except:",str(e)
        c,a=sv[-1],lag[-1] 
  return c,a 
#==============================================
# read inputs
CaMa_dir=sys.argv[1]
mapname=sys.argv[2]
inname=sys.argv[3]
out_dir=sys.argv[4]
para=int(sys.argv[5])
threshold=float(sys.argv[6])
#==============================================
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
with open(fname,"r") as f:
  lines=f.readlines()
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
print out_dir
#==============================================
# # spatial dependancy weight
# pathname0=out_dir+"/weightage"
# mk_dir(pathname0)
# print pathname0
# pathname1=out_dir+"/weightage/"+mapname+"_"+inname
# mk_dir(pathname1)
# print pathname1
#==============================================
# # gaussian weight for data assimilation
# pathname2=out_dir+"/gaussian_weight"
# mk_dir(pathname2)
# print pathname2
thresname="%02d"%(int(threshold*100))
pathname3=out_dir+"/gaussian_weight/"+mapname+"_"+inname+"_"+thresname
# mk_dir(pathname3)
print pathname3
#==============================================
fname=out_dir+"/semivar/"+mapname+"_"+inname+"/lonlat_list.txt"
with open(fname,"r") as f:
  lines=f.readlines()
#==============================================
# threshold for defining the local patch boundries
#threshold=0.6
# samllest patch size along the river (upstream/downstream)
#baseline=6
baseline=11 # Revel_etal,. (2019) proves 21 x 21 patch works well in upstreams 
#print lines[0]
#for line in lines[1::]:
#==============================================
def mk_wgt(line):
    #print line
    line  = filter(None, re.split(" ",line))
    # print line[0], line[1], line[2], line[3]
    lon   = int(line[0])
    lat   = int(line[1])
    up    = int(line[2])
    dn    = int(line[3])
    #---
    if up + dn == 0:
        return 0
    #--
    wgt=np.zeros([ny,nx],np.float32)
    Gwt=np.zeros([ny,nx],np.float32)
    #--
    # put 1 in self pixel (target pixel should be 1)
    wgt[lat-1,lon-1]=1.0
    Gwt[lat-1,lon-1]=1.0
    #---
    if dn>0:
        # print "downstream"
        #--
        lix=[]
        liy=[]
        ldis=[]
        lgamma=[]
        lstd=[]
        #--
        fname="%s/semivar/%s_%s/%04d%04d/dn%05d.svg"%(out_dir,mapname,inname,lon,lat,0)
        try:
            f=open(fname,"r")
            lines=f.readlines()
            f.close()
        except:
            print "no file", fname
            return 0
        #--
        for line in lines[1::]:
            line  = filter(None, re.split(" ",line))
            #print line
            ix    = int(line[0])
            iy    = int(line[1])
            dis   = float(line[2])
            gamma = float(line[3])
            std   = float(line[4])
            #print dis,gamma
            #--
            lix.append(ix)
            liy.append(iy)
            ldis.append(dis)
            lgamma.append(gamma)
            lstd.append(std)
        # find auto correlation length
        c,a=mk_svg(lgamma,ldis)
        #print c, a
        # cal weigtage
        lwgt=1.0 - np.array(lgamma)/(c+1.0e-20)
        #--
        if sum((lwgt>=threshold)*1)<baseline: #6
            a1=1.0e-20
            if len(ldis)>1:
                a1=(ldis[-1]/float(len(ldis)-1))*(float(baseline) - 1.)
                sigma=a1/math.sqrt(2.0*abs(math.log(threshold))) # calculate sigma for r=a1 and w-0.6 using w=exp(-r^2 /2*sigma^2)
                #define sigma => r = 2 * sqrt(10/3) * sigma
                sigma1=math.sqrt(3.0/10.0)*a1/2.0
                #---
                for i in range(0,len(ldis)-1):
                    x=lix[i]-1
                    y=liy[i]-1
                    #wgt[y,x]=max(gaussian(ldis[i],a),wgt[y,x])
                    wgt[y,x]=max(math.exp(-ldis[i]**2/(2.0*sigma**2)),wgt[y,x])
                    # Gaussian Weight
                    Gwt[y,x]=max(math.exp(-ldis[i]**2/(2.0*sigma1**2)),Gwt[y,x])
            else: # river mouth, inland termination
                if nextx[lat-1,lon-1] == -9 : # river mouth
                    x=lix[0]-1
                    y=liy[0]-1
                    wgt[y,x]=1.0
                    # Gaussian Weight
                    Gwt[y,x]=1.0
                #--
                if nextx[lat-1,lon-1] == -10 : # inland termination
                    x=lix[0]-1
                    y=liy[0]-1
                    wgt[y,x]=1.0
                    # Gaussian Weight
                    Gwt[y,x]=1.0
                # -9999
                if nextx[lat-1,lon-1] == -9999 :
                    print "-9999"
                    #return 0
        else:
            if sum((lwgt>=threshold)*1) >= len(ldis):
                a1=ldis[-1]
            else:
                a1=ldis[int(sum((lwgt>=threshold)*1))]
            #define sigma
            sigma1=math.sqrt(3.0/10.0)*a1/2.0
            for i in range(0,len(ldis)-1):
                x=lix[i]-1
                y=liy[i]-1
                #wgt[y,x]=max(gaussian(ldis[i],a),wgt[y,x])
                wgt[y,x]=max(1.0-(lgamma[i]/(c+1.0e-20)),wgt[y,x])
                # Gaussian Weight
                Gwt[y,x]=max(math.exp(-ldis[i]**2/(2.0*sigma1**2)),Gwt[y,x])
    #--
    #print "upstream"
    #plt.clf()
    if up > 0:
        # print "upstream"
        for iup in np.arange(1,up+1):
            #print "upstream",iup
            #--
            lix=[]
            liy=[]
            ldis=[]
            lgamma=[]
            lstd=[]
            #--
            fname="%s/semivar/%s_%s/%04d%04d/up%05d.svg"%(out_dir,mapname,inname,lon,lat,iup)
            try:
                f=open(fname,"r")
                lines=f.readlines()
                f.close()
            except:
                print "no file", fname
                #continue
                return 0
            #--
            for line in lines[1::]:
                line  = filter(None, re.split(" ",line))
                ix    = int(line[0])
                iy    = int(line[1])
                dis   = float(line[2])
                gamma = float(line[3])
                std   = float(line[4])
                #print dis,gamma
                #--
                lix.append(ix)
                liy.append(iy)
                ldis.append(dis)
                lgamma.append(gamma)
                lstd.append(std)
            #--
            c,a=mk_svg(lgamma,ldis)
            #print c, a
            # cal weigtage
            lwgt=1.0 - np.array(lgamma)/(c+1.0e-20)
            #--
            if sum((lwgt>=threshold)*1)<baseline:
                a1=1.0e-20
                if len(ldis)>1:
                    a1=(ldis[-1]/float(len(ldis)-1))*(float(baseline) - 1.) # calculate sigma for r=a1 and w-0.6 using w=exp(-r^2 /2*sigma^2)
                    #define sigma => r = 2 * sqrt(10/3) * sigma
                    sigma=a1/math.sqrt(2.0*abs(math.log(threshold)))
                    #define sigma
                    sigma1=math.sqrt(3.0/10.0)*a1/2.0
                    #---
                    for i in range(0,len(ldis)-1):
                        x=lix[i]-1
                        y=liy[i]-1
                        #wgt[y,x]=max(gaussian(ldis[i],a),wgt[y,x])
                        wgt[y,x]=max(math.exp(-ldis[i]**2/(2.0*sigma**2)),wgt[y,x])
                        # Gaussian Weight
                        Gwt[y,x]=max(math.exp(-ldis[i]**2/(2.0*sigma1**2)),Gwt[y,x])
                else: # only that grid
                        x=lix[0]-1
                        y=liy[0]-1
                        wgt[y,x]=1.0
                        # Gaussian Weight
                        Gwt[y,x]=1.0
            else:
                if sum((lwgt>=threshold)*1) >= len(ldis):
                    a1=ldis[-1]
                else:
                    a1=ldis[int(sum((lwgt>=threshold)*1))]
                #define sigma
                sigma1=math.sqrt(3.0/10.0)*a1/2.0
                for i in range(0,len(ldis)-1):
                    x=lix[i]-1
                    y=liy[i]-1
                    #wgt[y,x]=max(gaussian(ldis[i],a),wgt[y,x])
                    wgt[y,x]=max(1.0-(lgamma[i]/(c+1.0e-20)),wgt[y,x])
                    # Gaussian Weight 
                    Gwt[y,x]=max(math.exp(-ldis[i]**2/(2.0*sigma1**2)),Gwt[y,x])
    # #-----
    # oname=pathname1+"/%04d%04d.bin"%(lon,lat)
    # wgt.tofile(oname)
    # print oname
    #-----
    oname=pathname3+"/%04d%04d.bin"%(lon,lat)
    Gwt.tofile(oname)
    print oname
    return 0
##############################################
if para>0:
    print "do it parallel"
    os.system("export OMP_NUM_THREADS=%d"%(para))
    p=Pool(para)
    p.map(mk_wgt,lines[1::])
    p.terminate()
else:
    print "do it linear"
    map(mk_wgt,lines[1::])