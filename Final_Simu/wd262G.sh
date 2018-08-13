#!/bin/bash
#SBATCH --job-name=yiming

module load Apps/R
module load Tools/SimpleQueue ## for grace
sq_path=/SAY/archive/hz27-CC0937-Biostatistics-A/wd262/utmost_revision/sq_files
partition=day
for k in {8251..8500}
do
cd ${sq_path}/${k}
sqCreateScript -q ${partition} -N ${k} -n 60 -m 10000 -w 24:00:00 sq.sh > sq.sh.pbs
sbatch sq.sh.pbs
done