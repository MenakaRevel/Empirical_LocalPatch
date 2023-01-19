#!/bin/sh

#### initial setting ========================================
## project name
TAG=glb_15min

OUT_DAMPARAM="./${TAG}/damparam_${TAG}.csv"

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

mkdir -p ./${TAG}

echo "============================"
echo "DAM allocation project: $TAG"
echo "============================"

echo ""
echo "#####  Part 1: allocate GRanD on CaMa map"
echo ""
echo "@@@ src/get_rivinfo_glb "

./src/get_rivinfo_glb
mv ./damloc_tmp.txt ./$TAG/

# temporal dam allocation file : $TAG/damloc_tmp.txt

echo ""
echo "@@@ ./src/modify_damloc.py "
python ./src/modify_damloc.py $TAG $MINUPAREA

# mnodified dam allocation file : $TAG/damloc_modified.csv