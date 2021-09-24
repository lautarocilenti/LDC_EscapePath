#!/bin/bash
#SBATCH --job-name=Lautaro
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12




#SBATCH --mail-type=end
#SBATCH --mail-user=lcilenti@umd.edu

#### load and unload modules you may need
#module load cuda
#module load eog
#module load singularity
#module load tensorflow/1.10.1-gpu
#module load matlab/R2018a
module load matlab/R2019a
module list

SECONDS=0;
for i in {1..1}
do
export TMPDIR=$(pwd)/matlab_tmp # or select another path
mkdir -p $TMPDIR
matlab -nodisplay -nosplash -r "run('Main.m')"
done
echo "Seconds:";
echo $SECONDS;

echo "Finished with job $SLURM_JOBID"
