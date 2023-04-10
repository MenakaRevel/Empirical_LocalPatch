#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N mk_fig

#===========================
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

NCPUS0=20
export OMP_NUM_THREADS=$NCPUS0

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Empirical_LocalPatch/img"

# python mean_dis_map.py &

# python mean_wse_map.py &

# python map_weight.py

thresname="1000KM" 
damrep=1 
NCPUS=20
# python local_patch1.py $thresname $damrep $NCPUS

thresname="1000KM" 
damrep=0
NCPUS=20
# python local_patch1.py $thresname $damrep $NCPUS

thresname="60" 
damrep=1
NCPUS=20
python local_patch1.py $thresname $damrep $NCPUS

thresname="60" 
damrep=0
NCPUS=20
# python local_patch1.py $thresname $damrep $NCPUS

wait

conda deactivate