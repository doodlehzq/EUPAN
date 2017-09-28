s/LSF/SLURM/
s/lsf/slurm/

/BSUB -J/ {i\	    print JOB "\\#!/bin/bash\\n";
}
s/BSUB -J /SBATCH --job-name=/
s/BSUB -q /SBATCH -p /
s/BSUB -o /SBATCH --output=/
s/BSUB -e /SBATCH --error=/
s/BSUB -n /SBATCH -n /
s/BSUB -R .*ptile=\(.*\)]\\"/SBATCH --ntasks-per-node=\1/

s/bsub </sbatch /