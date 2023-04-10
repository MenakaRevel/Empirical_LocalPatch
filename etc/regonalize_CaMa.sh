#!/bin/sh
#====================
# Use externally simulated CaMa-Flood results to create Emperical Local Patches
# Convert global data to regional netCDF in same map resolution
# Better for high resoultion CaMa maps [e.g., 06min]
# From here one can start with s03-remove_trend.sh
# Menaka@IIS
# 2023/01/18
#====================
#*** PBS setting when needed
#PBS -q E20
#PBS -l select=1:ncpus=20:mem=100gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N regonalize
#===========================
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Empirical_LocalPatch/etc"

#================================================
# OpenMP Thread number
NCPUS=40
export OMP_NUM_THREADS=$NCPUS

# link params.py
rm -r params.py
ln -sf ../params.py params.py

# input settings
syear=2000 #`python -c "import params; print (params.starttime()[0])"`
smonth=`python -c "import params; print (params.starttime()[1])"`
sdate=`python -c "import params; print (params.starttime()[2])"`
eyear=2000 #`python -c "import params; print (params.endtime()[0])"`
emonth=`python -c "import params; print (params.endtime()[1])"`
edate=`python -c "import params; print (params.endtime()[2])"`
echo $syear" to "$eyear
CAMADIR=`python -c "import params; print (params.CaMa_dir())"`
outdir=`python -c "import params; print (params.out_dir())"`
# cpunums=`python -c "import params; print (params.cpu_nums())"`
cpunums=$NCPUS
mapname=`python -c "import params; print (params.map_name())"`
inputname=`python -c "import params; print (params.input_name())"`
N=`python src/calc_days.py $syear $smonth $sdate $eyear $emonth $edate`
threshold=`python -c "import params; print (params.threshold())"`

patch=100

# global map name
glbmapname="glb_06min"

# # link simulated files
# cd "../CaMa_out/${mapname}_${inputname}"
# pwd
# ln -sf "/home/yamadai/data/CaMa_v400_simulations/VIC_BC_3h_06min/* ."
# cd ../../etc

#======
# create folder name in CaMa_out
mkdir -p "../CaMa_out/${mapname}_${inputname}"

# water surface elevation
varname="sfcelv"
#=================================================
echo ./src/bin2nc_reg $N $syear $eyear $varname $mapname $glbmapname $inputname $CAMADIR $outdir &
time ./src/bin2nc_reg $N $syear $eyear $varname $mapname $glbmapname $inputname $CAMADIR $outdir &

# discharge
varname="outflw"
#=================================================
echo ./src/bin2nc_reg $N $syear $eyear $varname $mapname $glbmapname $inputname $CAMADIR $outdir &
time ./src/bin2nc_reg $N $syear $eyear $varname $mapname $glbmapname $inputname $CAMADIR $outdir &

wait