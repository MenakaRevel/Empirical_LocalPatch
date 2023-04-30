#!/bin/sh

# This code is modefied from CaMa-Flood v4.07
# written by Risa Hanazaki & Dai Yamazaki @UTokyo

## input data
#- GRanD_v1_1_inp.csv
#[damid, damname, lon, lat, total storage capacity, drainage area]

#### initial setting ========================================
## project name
TAG="glb_15min"
# TAG="glb_06min"

GRANDF="/home/yamadai/work/data/Dam+Lake/GRanD/inp/GRanD_v1_1_inp.csv"

OUT_TMP="../dat/damloc_tmp.txt"

# OUT_DAMLOC="../dat/damloc_modified.csv"
OUT_DAMLOC="../dat/damloc_${TAG}.txt"

MAPDIR="/cluster/data6/menaka/CaMa-Flood_v4/map/${TAG}"

## minimum drainage area [km2]
MINUPAREA=1000

## Naturalized simulation (natsim) parameters
#Sample test1: test1-glb_15min
SYEAR=2000   ## Start year of natsim
EYEAR=2001   ## End   year of natsim
DT=86400     ## time step of outflw.bin data (sec)

#Sample test4: test4-e2o_ecmwf-glb_15min
#SYEAR=1980   ## Start year of natsim
#EYEAR=2014   ## End   year of natsim
#DT=86400     ## time step of outflw.bin data (sec)

#===========================================================

# mkdir -p ./${TAG}

echo "============================"
echo "DAM allocation project: $TAG"
echo "============================"

echo ""
echo "#####  allocate GRanD on CaMa map"

# ./t01-calc_damloc.sh   $TAG $MINUPAREA

# output dam allocation file : $TAG/damloc_modified.csv

##==========================

echo ""
echo "@@@ src/get_rivinfo_glb "
echo ./src/get_rivinfo_glb ${GRANDF} ${OUT_TMP} ${MAPDIR}
./src/get_rivinfo_glb ${GRANDF} ${OUT_TMP} ${MAPDIR}


echo ""
echo "@@@ ./src/modify_damloc.py "
python ./src/modify_damloc.py $TAG $MINUPAREA $MAPDIR $OUT_TMP $OUT_DAMLOC
