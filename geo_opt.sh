#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --error=err-%j.log
#SBATCH --output=out-%j.log
#SBATCH --time=3-22:00:00

 
HDIR=$(pwd)

export WDIR=$TMPDIR/$SLURM_JOB_NAME

srun mkdir -p $WDIR

#cpc $WDIR
cp $HDIR/* $WDIR
cd $WDIR


export TMP=$TMPDIR/$SLURM_JOBID
# do not remove /$SLURM_JOBID
 
srun mkdir -p $TMP
# -p option make the $TMPDIR directory on all nodes
  
#setup turbomole:
export PARA_ARCH=SMP
export PARNODES=20
export TURBODIR=/opt/software/user-soft/turbomole/76
export PATH=$TURBODIR/bin/`sysname`:$PATH
export PATH=$TURBODIR/scripts:$PATH
export TURBOTMPDIR=$TMP
  
# TURBOMOLE command(s):
jobex -gcart 4 -c 1000 > jobex.out
 
rm -rf auxboi auxbai ftraux fock oldfock errvec
rm -rf dens ddens mdens core twoint*
rm -rf CC_* CC1DEN* CCJCOUINTE CCTR* CCVPQ_ISQR CCB* CCG* CCEQB* CCY* CCI* CCFOCK_LMO
rm -rf RIR12_*
rm -rf half* gamma* moint* alles.cc syminfo oneint
rm -rf statistics.*

rm -rf $TMP*

mkdir $HDIR/BACK
cp -r $WDIR/* $HDIR/BACK/.
