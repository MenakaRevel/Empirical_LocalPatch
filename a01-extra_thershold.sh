#!/bin/sh
#====================
# Calcualte the experimental semivarinces along each river stem
# Menaka@IIS
# 2020/06/01
#====================
#*** PBS setting when needed
#PBS -q F40
#PBS -l select=1:ncpus=40:mem=40gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N e_threshold

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
# threshold=`python -c "import params; print (params.threshold())"`
threshold=0.92

patch=100

threshname=$(echo $threshold 100 | awk '{printf "%2d\n",$1*$2}')

#=================================================
# Gaussian weight for given thershold

# make dir local patch
mkdir -p "./gaussian_weight/${mapname}_${inputname}_${threshname}"

# Make Gaussian weight
python src/gaussian_weight.py $CAMADIR $mapname $inputname $outdir $cpunums $threshold

#=================================================
# Full local patch
# make dir local patch
mkdir -p "./local_patch/${mapname}_${inputname}_${threshname}"

# Write local patch parameters
./src/lpara $N $syear $eyear $mapname $inputname $CAMADIR $outdir $threshold $patch

#=================================================
# Mainstream local patch
# make dir local patch
mkdir -p "./local_patchMS/${mapname}_${inputname}_${threshname}"

# Write local patch parameters
./src/lparaMS $N $syear $eyear $mapname $inputname $CAMADIR $outdir $threshold

wait

conda deactivate