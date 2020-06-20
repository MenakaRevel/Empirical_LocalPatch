#!/opt/local/bin/python
import numpy as np
import os
#! /usr/bin/python
import re
import numpy as np
"""give the edges of south north west east corrdinates for mapping for major rivers
Menaka@IIS 2019/11/01"""
def latlon_river(rivername):
     global lllat, urlat, lllon, urlon
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
     if rivername=="AMAZONAS":
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
         lllat = 30.
         urlat = 50.
         lllon = -110.
         urlon = -85.
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
     return lllat, urlat, lllon, urlon
#-------------------------------------
