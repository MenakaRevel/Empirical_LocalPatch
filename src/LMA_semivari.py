from pylab import *
import numpy as np
from scipy.spatial.distance import pdist, squareform
import math
from scipy import *
from scipy.optimize import leastsq
import re
import sys, os

# Levenberg-Marquardt method
# **************************
def initial(x,y):
    C0i=max(0.0,y[0]-((x[0]/(x[1]-x[2]))*(y[1]-y[0])))
    #--
    ai=x[-1]/2.0
    #--
    Ci=(sum(y[-3::])/3.0)-C0i
    #--
    return C0i,ai,Ci
#----------------
def residuals(p,model_name, y, x,fit):
    model={'spherical':spherical, 'cubic':cubic, 'exponential':exponential, 'gaussian':gaussian, 'pentaspherical':pentaspherical, 'sineholeeffect':sineholeeffect,'power':power}
    err = model[model_name](x,p,fit) - y
    return err
#------------------------
def mrqfit(model_name,x,y,fit=[1,1,1]):
    C0i, ai, Ci=initial(x,y)
    p0=array([C0i, ai, Ci])
    #--
    if sum(fit)==3: 
      plsq = leastsq(residuals, p0, args=(model_name,y, x,fit),Dfun=Jacobian, maxfev=2000)
#      plsq = leastsq(residuals, p0, args=(model_name,y, x),maxfev=2000)
    elif fit[0]==0:
      plsq = leastsq(residuals, p0[1:], args=(model_name,y, x,fit),Dfun=Jacobian, maxfev=2000)

    return plsq
#------------------------
def fit_para(p,fit):
    p0=np.array([0.,0.,0.])
    if sum(fit) < 3:
      p0[0]=0.
      p0[1]=p[0]
      p0[2]=p[1]
    else:
      p0=p
    return p0
#-------------------------
def AIC_val(model,plsq,n,x,y,fit=[1,1,1]):
    p=float(sum(fit))
    resid = sum(np.sqrt((y-pval(model,x,plsq,fit))**2))
    #n=float(n) 
    n=float(len(y))
    #--
    AIC=n*math.log(resid/n)+2*p
    return AIC
#--------------------------
def pval(model_name,x,p,fit=[1,1,1]):
    model={'spherical':spherical, 'cubic':cubic, 'exponential':exponential, 'gaussian':gaussian, 'pentaspherical':pentaspherical, 'sineholeeffect':sineholeeffect,'power':power}
    return model[model_name](x,p,fit)
#--------------------------
def Jacobian(p,model_name, y, x,fit=[1,1,1]):
    Jmodel={'spherical':Jspherical, 'cubic':Jcubic, 'exponential':Jexponential, 'gaussian':Jgaussian, 'pentaspherical':Jpentaspherical, 'sineholeeffect':Jsineholeeffect,'power':Jpower}
    return Jmodel[model_name](x,p,fit)
#--------------------------
# Semi-Variogram Models
#*************************
def H(n):
    if n==0:
      h=1
    else:
      h=0
    return h
#---------------------------
def spherical(x, p,fit=[1,1,1]):
    p=fit_para(p,fit)
    yfit=[]
    for i in x:
      if i<p[1]:
        yfit.append( (p[0]*(1.-H(i)) + p[2]*((3./2.)*(i/p[1])-(1./2.)*(i/p[1])**3.)) )
      elif i>=p[1]:
        yfit.append( (p[0] + p[2]) )
    return np.array(yfit)
#---------------
def cubic(x, p,fit=[1,1,1]):
    p=fit_para(p,fit)
    yfit=[]
    for i in x:
      if i<p[1]:
        yfit.append( (p[0]*(1.-H(i)) + p[2]*((7.0)*(i/p[1])**2.0-(35./4.)*(i/p[1])**3.+(7./2.)*(i/p[1])**5.-(3./4.)*(i/p[1])**7.)) )
      elif i>=p[1]:
        yfit.append( (p[0] + p[2]) )
    return np.array(yfit)
#------------------------
def exponential(x, p,fit=[1,1,1]):
    p=fit_para(p,fit)
    yfit=[]
    for i in x:
      if i>=0:
        yfit.append( (p[0]*(1.-H(i)) + p[2]*(1.0 - math.exp(-(3.0*i)/p[1]))))
    return np.array(yfit)
#------------------------
def power(x, p,fit=[1,1,1]):
    p=fit_para(p,fit)
    yfit=[]
    for i in x:
      if i>=0:
        yfit.append( (p[0] * (1.-H(i)) + p[2] * (i**p[1])))
    return np.array(yfit)
#------------------------
def gaussian(x, p,fit=[1,1,1]):
    p=fit_para(p,fit)
    yfit=[]
    for i in x:
      if i>=0:
        yfit.append( (p[0] * (1.-H(i)) + p[2]*(1.0 - math.exp(-3.0*(i/p[1])**2.0))))
    return np.array(yfit)
#------------------------
def pentaspherical(x, p,fit=[1,1,1]):
    p=fit_para(p,fit)
    yfit=[]
    for i in x:
      if i<p[1]:
        yfit.append( (p[0]*(1.-H(i)) + p[2]*((15.0/8.0)*(i/p[1])-(5./4.)*(i/p[1])**3.+(3./8.)*(i/p[1])**5.)) )
      elif i>=p[1]:
        yfit.append( (p[0] + p[2]) )
    return np.array(yfit)
#------------------------
def sineholeeffect(x, p,fit=[1,1,1]):
    p=fit_para(p,fit)
    yfit=[]
    for i in x:
      if i>=0:
        yfit.append( (p[0] * (1.-H(i)) + p[2] * (1.0 - math.sin(4.4934*(i/p[1]))/(4.4934*(i/p[1])))) )
    return np.array(yfit)
#---------------------------
def Jspherical(x,p,fit=[1,1,1]):
    p=fit_para(p,fit)
    dy=np.zeros([len(x),3],np.float32)
    for n,i in enumerate(x):
      if i<p[1]:
        dy[n,0]= (1.-H(i)) 
        dy[n,1]= (p[2]*(-(3./2.)*(i/p[1])+3.*(1./2.)*(i/p[1])**3.))/p[1]
        dy[n,2]= ((3./2.)*(i/p[1])-(1./2.)*(i/p[1])**3.)
      elif i>=p[1]:
        dy[n,0]= 1.0
        dy[n,1]= 0.0
        dy[n,2]= 1.0
    if sum(fit)<3:
      return np.array(dy[:,1:])
    else:
      return np.array(dy)
#---------------
def Jcubic(x,p,fit=[1,1,1]):
    p=fit_para(p,fit)
    dy=np.zeros([len(x),3],np.float32)
    for n,i in enumerate(x):
      if i<p[1]:
        dy[n,0]= (1.-H(i))
        dy[n,1]= (p[2]*(-2.*(7.)*(i/p[1])**2.0+3.*(35./4.)*(i/p[1])**3.-5.*(7./2.)*(i/p[1])**5.+7.*(3./4.)*(i/p[1])**7.))/p[1]
        dy[n,2]= ((7.0)*(i/p[1])**2.0-(35./4.)*(i/p[1])**3.+(7./2.)*(i/p[1])**5.-(3./4.)*(i/p[1])**7.)
      elif i>=p[1]:
        dy[n,0]= 1.0
        dy[n,1]= 0.0
        dy[n,2]= 1.0
    if sum(fit)<3:
      return np.array(dy[:,1:])
    else:
      return np.array(dy)

#------------------------
def Jexponential(x,p,fit=[1,1,1]):
    p=fit_para(p,fit)
    dy=np.zeros([len(x),3],np.float32)
    for n,i in enumerate(x):
      if i>=0:
        dy[n,0]= (1.-H(i))
        dy[n,1]= (p[2]*(- math.exp(-(3.0*i)/p[1])*(i/p[1])))/p[1]
        dy[n,2]= (1.0 - math.exp(-(3.0*i)/p[1])) 
    if sum(fit)<3:
      return np.array(dy[:,1:])
    else:
      return np.array(dy)

#------------------------
def Jpower(x,p,fit=[1,1,1]):
    p=fit_para(p,fit)
    dy=np.zeros([len(x),3],np.float32)
    for n,i in enumerate(x):
      if i>=0:
        dy[n,0]= (1.-H(i))
        dy[n,1]= (p[2]*(i**p[1])*math.log(i))
        dy[n,2]= (i**p[1]) 
    if sum(fit)<3:
      return np.array(dy[:,1:])
    else:
      return np.array(dy)
#------------------------
def Jgaussian(x,p,fit=[1,1,1]):
    p=fit_para(p,fit)
    dy=np.zeros([len(x),3],np.float32)
    for n,i in enumerate(x):
      if i>=0:
        dy[n,0]= (1.-H(i))
        dy[n,1]= (p[2]*(math.exp(-3.0*(i/p[1])**2.0))*(-3.*2.*(i/p[1])**2.0))/p[1]
        dy[n,2]= (1.0 - math.exp(-3.0*(i/p[1])**2.0))
    if sum(fit)<3:
      return np.array(dy[:,1:])
    else:
      return np.array(dy)
#------------------------
def Jpentaspherical(x,p,fit=[1,1,1]):
    p=fit_para(p,fit)
    dy=np.zeros([len(x),3],np.float32)
    for n,i in enumerate(x):
      if i<p[1]:
        dy[n,0]= (1.-H(i))
        dy[n,1]= (p[2]*(-(15.0/8.0)*(i/p[1])+3.*(5./4.)*(i/p[1])**3.-5.*(3./8.)*(i/p[1])**5.))/p[1]
        dy[n,2]= ((15.0/8.0)*(i/p[1])-(5./4.)*(i/p[1])**3.+(3./8.)*(i/p[1])**5.)
      elif i>=p[1]:
        dy[n,0]= 1.0
        dy[n,1]= 0.0
        dy[n,2]= 1.0
    if sum(fit)<3:
      return np.array(dy[:,1:])
    else:
      return np.array(dy)
#------------------------
def Jsineholeeffect(x,p,fit=[1,1,1]):
    p=fit_para(p,fit)
    dy=np.zeros([len(x),3],np.float32)
    for n,i in enumerate(x):
      if i>=0:
        dy[n,0]= (1.-H(i))   
        dy[n,1]= (p[2]*math.cos(i)*(math.sin(4.4934*(i/p[1]))/(4.4934*(i/p[1]))))/p[1] 
        dy[n,2]= (1.0 - math.sin(4.4934*(i/p[1]))/(4.4934*(i/p[1])))
    if sum(fit)<3:
      return np.array(dy[:,1:])
    else:
      return np.array(dy)
#------------------------

