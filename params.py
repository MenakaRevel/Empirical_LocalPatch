#!/opt/local/bin/python
# -*- coding: utf-8 -*-
# import numpy as np
########################
# Creating Empirical Local Patch
# By Menaka Revel@IIS,U-Tokyo
# 2020/05/25
# Revel et al,. (2019,2021,2022)
########################
#
# parameters list
#
########################

def timestep():
     return 86400 # outer timestep in seconds

def starttime():
     return [1979,1,1] # start date: [year,month,date]

def endtime():
     return [2019,12,31] # end date: [year,month,date]
                      # *note: this date is not included

def CaMa_dir():
    #return "/cluster/data6/menaka/CaMa-Flood_v395b_20191030"
    #return "/cluster/data6/menaka/CaMa-Flood_v396_20191225"
    # return "/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
    return "/cluster/data6/menaka/CaMa-Flood_v4"
    # directory of CaMa-Flood
    # indicate the directory of ./map or ./src and other folders

def out_dir():
    return "/cluster/data6/menaka/Empirical_LocalPatch"
    # return "/work/a06/menaka/Empirical_LocalPatch"
    #return "/media/menaka/HDJA-UT/covariance"

def map_name():
    # return "amz_06min" #Amazon river 6 min
    # return "congo_06min" #Congo river 6min
    # return "glb_06min" # global 6 min
    return "glb_15min" # global 15 min

def input_name():
    return "S14FD"
    # return "VIC_BC"

def spinup_mode():
     return 0
     # 0: no spinup simulation 
     # 1: do spinup simulation
     ### if initial restart file is ready, spinup simulation is no need

def spinup_end_year():
    return 1979

def spinup_end_month():
    return 12

def spinup_end_date():
    return 31

def patch_start():
    return 1979,1,1

def patch_end():
    return 2019,12,31

def threshold():
    return 0.6000000
    # return 0.4000000
    # return 0.2000000

def dam_rep():
    return 1 # represent dams in local patch
    # return 0 # represent dams in local patch

def para_nums():
    return 6
    # setting number of parallels to run CaMa-Flood Model
    # defualt is 6, but may change depending on your system

def slack_notification():
    return 0
    # setting for validating slack notification
    # 1 for valid and 0 for invalid
    # 0 is a default if you are not familiar with slack
    # if you turn it to 1, you need to edit sendslack.py
    # for more information refer https://api.slack.com/incoming-webhooks

def cpu_nums():
    return 20
    # number of cpus used 

# def qsub():
#     return 0
#     # 1: submit the job via qsub , resource allocation should be done
#     # 0: direct execution in forntnode

# def qoption():
#     return "-q E40 -l select=1:ncpus=40:mem=10gb -d "
#     # if qsub = 1 the theses options will be used

def version():
    # return "v3.0.0 (updated 2020-05-25): CaMa-Flood v396"
    return "v3.0.1 (updated 2022-11-23): CaMa-Flood v4, dams included"
    # version  396 merit DEM
    # differnt maps [glb_15min,glb_06min, etc]
