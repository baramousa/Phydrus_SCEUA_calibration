#!/bin/bash -l
#SBATCH --job-name=test_serial
#SBATCH --ntasks=14
#SBATCH --mem-per-cpu=3G
#SBATCH --time=200:00:00
#SBATCH --constraint=[skylake|haswell]
#SBATCH --mail-user=albara@uni-hannover.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output test_serial-job_%j.out
#SBATCH --error test_serial-job_%j.err

module load GCCcore/.9.3.0 Python/3.8.2
module load GCC/9.3.0 OpenMPI/4.0.3  numpy/1.18.5-Python-3.8.2
module load mpi4py/3.0.3-Python-3.8.2

plot="ir_lablab_med"
home_plot="$HOME/spotpy_hyd/"$plot
bigwork="$BIGWORK/spotpy_hyd"
bigwork_plot=$bigwork$"/"$plot
bigwork_source=$bigwork_plot$"/source"

cp -r $home_plot $bigwork
cd $bigwork_source
make
sleep 10
cp hydrus $bigwork_plot
cd $bigwork_plot
srun python execute.py
sleep 30
cp -r $bigwork_plot $HOME/spotpy_hyd/results
cp -r $bigwork$"/"$plot$"opt/re_17_test" $HOME/spotpy_hyd/results/$plot