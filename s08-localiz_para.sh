#!/bin/sh
#====================
# Write localization parameters to easy acess text files
# Menaka@IIS
# 2020/06/01
#====================
#*** PBS setting when needed
#PBS -q F10
#PBS -l select=1:ncpus=10:mem=40gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N lpara
#========
cd $PBS_O_WORKDIR
#================================================
# OpenMP Thread number
export OMP_NUM_THREADS=10

# input settings
syear=`python -c "import params; print (params.starttime()[0])"`
smonth=`python -c "import params; print (params.starttime()[1])"`
sdate=`python -c "import params; print (params.starttime()[2])"`
eyear=`python -c "import params; print (params.endtime()[0])"`
emonth=`python -c "import params; print (params.endtime()[1])"`
edate=`python -c "import params; print (params.endtime()[2])"`
echo $syear" to "$eyear
CAMADIR=`python -c "import params; print (params.CaMa_dir())"`
# CAMADIR="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
outdir=`python -c "import params; print (params.out_dir())"`
cpunums=`python -c "import params; print (params.cpu_nums())"`
mapname=`python -c "import params; print (params.map_name())"`
# mapname="amz_06min" #
inputname=`python -c "import params; print (params.input_name())"`
N=`python src/calc_days.py $syear $smonth $sdate $eyear $emonth $edate`
# threshold=`python -c "import params; print (params.threshold())"`
threshold=0.60
patch=100

threshname=$(echo $threshold 100 | awk '{printf "%2d\n",$1*$2}')

# make dir local patch
mkdir -p "./local_patch/${mapname}_${inputname}_${threshname}"

#=================================================
varname="weightage"
./src/lpara $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir $threshold $patch
