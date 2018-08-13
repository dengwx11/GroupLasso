#!/bin/bash
#SBATCH --partition=pi_zhao
#SBATCH --job-name=SNR10

module load Apps/R/3.3.2-generic
module load Tools/SimpleQueue ## for grace

Rscript final_simu.R -m_G 100