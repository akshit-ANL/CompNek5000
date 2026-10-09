#!/bin/bash -x

case="base"
jobename="det2d"
PROJ_ID="ML-GA"
QUEUE="compute"
nodes=20
procs=128
nn=$((nodes*procs))
time="24:00:00"

# sbatch
SFILE=s.bin
echo "#!/bin/bash -l" > $SFILE
echo " " >> $SFILE
echo "#PBS -A $PROJ_ID" >> $SFILE
echo "#PBS -N $jobename" >> $SFILE
echo "#PBS -q $QUEUE" >> $SFILE
echo "#PBS -l select=$nodes:mpiprocs=$procs" >> $SFILE
echo "#PBS -l walltime=$time" >> $SFILE
echo "#PBS -j eo" >> $SFILE
echo " " >> $SFILE
echo "cd \$PBS_O_WORKDIR" >> $SFILE
echo "echo $case > SESSION.NAME" >> $SFILE
echo "echo `pwd`/ >> SESSION.NAME" >> $SFILE
echo "module list" >> $SFILE
echo " " >> $SFILE
echo "mpiexec -n $nn ./nek5000 > $case.log 2>&1" >> $SFILE
qsub $SFILE



