ofname=`echo $1 | sed -e 's/\(.*\)LSF.pm/\1SLURM.pm/'`
sed -f lsf2slurm.sed $1 > $ofname