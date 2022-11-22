#!/bin/sh
#====================
# convert binary files to netCDF4 
# all the time steps will be in one file
# dimesnion (time, nx, ny)
# Menaka@IIS
# 2020/05/29
#====================
#*** PBS setting when needed
#PBS -q F40
#PBS -l select=1:ncpus=40:mem=100gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N EmpLP
#===========================
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
export OMP_NUM_THREADS=40

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

#*************************
# runs s01-CaMa_sim.sh before following steps to get the simulated variables
# ./s01-CaMa_sim.sh
#*************************
# # #=========================
# # # convert binary to netCDF
# # #=========================
# # # water surface elevation
# # varname="sfcelv"
# # #=================================================
# # ./src/bin2nc $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir

# # # discharge
# # varname="outflw"
# # #=================================================
# # ./src/bin2nc $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir

#=========================
# remove trend lines
#=========================
varname="sfcelv"
#=================================================
./src/remove_trend $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir

#=========================
# remove seasonality
#=========================
varname="rmdtrnd"
#=================================================
./src/remove_season $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir

#=========================
# standardized sfcelve
#=========================
varname="rmdsesn"
#=================================================
./src/standardize $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir

#=========================
# calculate experimental semi-varaince
#=========================
varname="standardized"
#===make directories for semivar
`python src/make_semivari.py $CAMADIR $mapname $outdir`
#=================================================
./src/semivariance $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir

#=========================
# calculate spatial auto-correlation weightage
#=========================
threshold=`python -c "import params; print (params.threshold())"`
threshname=$(echo $threshold 100 | awk '{printf "%2d\n",$1*$2}')
#=================================================
python src/weightage.py $CAMADIR $mapname $outdir $cpunums $threshold

#=========================
# write local patch to text files
#=========================
varname="weightage"

# make dir local patch
mkdir -p "local_patch/${mapname}_${inputname}_${threshname}"

#=================================================
./src/lpara $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir $threshold

#=========================
# write local patch [main stream] to text files
#=========================
varname="weightage"

# make dir local patch
mkdir -p "local_patchMS/${mapname}_${inputname}_${threshname}"
#=================================================
#./src/lparaMS $N $syear $eyear $varname $mapname $inputname $CAMADIR $outdir $threshold

wait

conda deactivate