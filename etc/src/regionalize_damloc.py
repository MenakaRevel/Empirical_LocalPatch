import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import re
#=======================
glbmap     = sys.argv[1]
regmap     = sys.argv[2]
inputfile  = sys.argv[3]
outputfile = sys.argv[4]
CaMa_dir   = sys.argv[5]
#=======================
# regional map
#==============================================
fname=CaMa_dir+"/map/"+regmap+"/params.txt"
with open(fname,"r") as f:
  lines=f.readlines()
#==============================================
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

# open global dam location file
dam_data = pd.read_csv(inputfile, sep='\s+', header=0)
print (dam_data.head())
for col in dam_data.columns:
  print (col)
# write text file 
with open(outputfile,"w") as fout:
    head="%6s%43s%10s%10s%8s%8s%12s%12s\n"%("Dam_ID","DamName","DamLon","DamLat","DamIX","DamIY","UpArea","Capacity")
    fout.write(head)
    for index, row in dam_data.iterrows():
        # print (row) #row["Dam_ID"],row['DamIX'])
        if row['DamIX']-offsetX < 0 or  row['DamIY']-offsetY < 0:
            continue      
        # damid   damname   damlon   damlat   ix   iy   upreal   uparea_cama   totalsto_mcm
        line="%06d%43s%10.2f%10.2f%8d%8d%12.2f%12.2f\n"%(row['Dam_ID'],row['DamName'],row['DamLon'],row['DamLat'],row['DamIX']-offsetX,row['DamIY']-offsetY,row['UpArea'],row['Capacity'])
        print(line)
        fout.write(line)