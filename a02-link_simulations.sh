#!/bin/sh
#====================
# Use externally simulated CaMa-Flood results to create Emperical Local Patches
# Menaka@IIS
# 2022/07/28
#====================
#*** PBS setting when needed
#PBS -q E20
#PBS -l select=1:ncpus=20:mem=100gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N e_link_sim

#===========================
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

#===========================
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Empirical_LocalPatch"

#================================================
# OpenMP Thread number
NCPUS=40
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
outdir=`python -c "import params; print (params.out_dir())"`
# cpunums=`python -c "import params; print (params.cpu_nums())"`
cpunums=$NCPUS
mapname=`python -c "import params; print (params.map_name())"`
inputname=`python -c "import params; print (params.input_name())"`
N=`python src/calc_days.py $syear $smonth $sdate $eyear $emonth $edate`
threshold=`python -c "import params; print (params.threshold())"`

patch=100

#======
# create folder name in CaMa_out
mkdir -p "./CaMa_out/${mapname}_${inputname}"

# link simulated files
# cd "./CaMa_out/${mapname}_${inputname}"
# ln -sf "/home/yamadai/data/CaMa_v400_simulations/VIC_BC_3h_06min/* ."
# cd ../../

# water surface elevation
varname="sfcelv"
#=================================================
./src/bin2nc $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir

# discharge
varname="outflw"
#=================================================
./src/bin2nc $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir