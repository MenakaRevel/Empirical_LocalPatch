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
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python
#========
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Empirical_LocalPatch"
#================================================
# OpenMP Thread number
NCPUS=10
export OMP_NUM_THREADS=$NCPUS

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
# represnt dams
damrep=0 #`python -c "import params; print (params.dam_rep())"`
# mapname="amz_06min" #
inputname=`python -c "import params; print (params.input_name())"`
N=`python src/calc_days.py $syear $smonth $sdate $eyear $emonth $edate`
threshold=`python -c "import params; print (params.threshold())"`
# threshold=0.60
patch=100
# distpatch=1 # distance based local patch
distpatch=0 # distance based local patch
threshname=$(echo $threshold 100 | awk '{printf "%2d\n",$1*$2}') # emperical local patch

# # local patch name
# if [ ${distpatch} -eq 1 ]; then
#     threshname="1000KM" # for distance based
# else
#     threshname=$(echo $threshold 100 | awk '{printf "%2d\n",$1*$2}') # emperical local patch
# fi

# make dir local patch
if [ ${damrep} -eq 1 ]; then
    mkdir -p "./local_patch/${mapname}_${inputname}_${threshname}_dam"
else
    mkdir -p "./local_patch/${mapname}_${inputname}_${threshname}"
fi

#=================================================
# Write local patch parameters
echo "./src/lpara $N $syear $eyear $mapname $inputname $CAMADIR $outdir $threshold $patch $damrep $NCPUS"
./src/lpara $N $syear $eyear $mapname $inputname $CAMADIR $outdir $threshold $patch $damrep $NCPUS

wait

conda deactivate