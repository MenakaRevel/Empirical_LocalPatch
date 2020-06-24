#!/bin/sh
#====================
# convert binary files to netCDF4 
# all the time steps will be in one file
# dimesnion (time, nx, ny)
# Menaka@IIS
# 2020/05/29
#====================
#*** PBS setting when needed
#PBS -q E20
#PBS -l select=1:ncpus=20:mem=100gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N bin2nc
#========
cd $PBS_O_WORKDIR
#================================================
# OpenMP Thread number
export OMP_NUM_THREADS=20

# input settings
syear=`python -c "import params; print (params.starttime()[0])"`
smonth=`python -c "import params; print (params.starttime()[1])"`
sdate=`python -c "import params; print (params.starttime()[2])"`
eyear=`python -c "import params; print (params.endtime()[0])"`
emonth=`python -c "import params; print (params.endtime()[1])"`
edate=`python -c "import params; print (params.endtime()[2])"`
echo $syear" to "$eyear
CAMADIR=`python -c "import params; print (params.CaMa_dir())"`
outdir=`python -c "import params; print (params.out_dir())"`
cpunums=`python -c "import params; print (params.cpu_nums())"`
mapname=`python -c "import params; print (params.map_name())"`
inputname=`python -c "import params; print (params.input_name())"`
N=`python src/calc_days.py $syear $smonth $sdate $eyear $emonth $edate`

# water surface elevation
varname="sfcelv"
#=================================================
./src/bin2nc $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir

# discharge
varname="outflw"
#=================================================
./src/bin2nc $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir
