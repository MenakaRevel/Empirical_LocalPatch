#!/bin/sh
#====================
# Calcualte the weights considering the distacne along each river stem [e.g., 500km]
# Menaka@IIS
# 2020/06/01
#====================
#*** PBS setting when needed
#PBS -q E40
#PBS -l select=1:ncpus=40:mem=40gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N dist_wgt
#========
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python
#========
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Empirical_LocalPatch/etc"
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
cpunums=$NCPUS #`python -c "import params; print (params.cpu_nums())"`
mapname=`python -c "import params; print (params.map_name())"`
inputname=`python -c "import params; print (params.input_name())"`
N=`python src/calc_days.py $syear $smonth $sdate $eyear $emonth $edate`
threshold=1000 #km distance based local patches ###`python -c "import params; print (params.threshold())"`
# represnt dams
damrep=1 #`python -c "import params; print (params.dam_rep())"`
#=================================================
python src/weightage_distance.py $CAMADIR $mapname $inputname $outdir $cpunums $threshold $damrep &

wait

conda deactivate