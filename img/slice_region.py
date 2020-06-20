from numpy import *
import numpy.ma as ma
#********************
#-----------------------------
def slice(data,lllat,urlat,lllon,urlon,res=1.0):
  if res == 1.0:  
    x1 = int((lllon + 0.5 - (-179.5))/1.0)
    x2 = int((urlon - 0.5 - (-179.5))/1.0)  
    y1 = int((lllat + 0.5 - (89.5))/1.0  )
    y2 = int((urlat - 0.5 - (89.5))/1.0  ) 
  elif res == 0.25:     
    x1 = int((lllon + 0.125 - (-179.875))/0.25)
    x2 = int((urlon - 0.125 - (-179.875))/0.25) 
    y1 = int(((89.875) - (urlat + 0.125))/0.25)
    y2 = int(((89.875) - (lllat - 0.125))/0.25)
  #--   
  data_region = data[y1:y2 + 1, x1:x2 + 1]  
  return data_region
#-----------------------------
def slice_aphro(data,lllat,urlat,lllon,urlon):
  x1 = (lllon + 0.125 - 60.125)/0.25
  x2 = (urlon - 0.125 - 60.125)/0.25 
  y1 = (lllat + 0.125 - (-14.875))/0.25
  y2 = (urlat - 0.125 - (-14.875))/0.25
  #--   
  data_region = data[y1:y2 + 1, x1:x2 + 1]  
  return data_region
#-----------------------------
def slice_GSMaP(data,lllat,urlat,lllon,urlon):
  x1 = (lllon + 0.005 - 0.005)/0.01
  x2 = (urlon - 0.005 - 0.005)/0.01 
  y1 = (lllat + 0.005 - (-59.995))/0.01
  y2 = (urlat - 0.005 - (-59.995))/0.01
  #--   
  data_region = data[y1:y2 + 1, x1:x2 + 1]  
  return data_region
#-----------------------------
def slice_SWOT(data,lllat,urlat,lllon,urlon,res=1.0):
  if res == 1.0:  
    x1 = (lllon + 0.5 - (-179.5))/1.0
    x2 = (urlon - 0.5 - (-179.5))/1.0  
    y1 = (lllat + 0.5 - (-79.5))/1.0
    y2 = (urlat - 0.5 - (-79.5))/1.0   
  elif res == 0.25:     
    x1 = (lllon + 0.125 - (-179.875))/0.25
    x2 = (urlon - 0.125 - (-179.875))/0.25 
    y1 = (lllat + 0.125 - (-79.875))/0.25
    y2 = (urlat - 0.125 - (-79.875))/0.25
  #--   
  data_region = data[y1:y2 + 1, x1:x2 + 1]  
  return data_region
