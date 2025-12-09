#! /bin/bash

for t in 0.2000 1.0000 
#for t in 5.0000
do
for h in 0.2000 1.0000
#for h in 1.0000
do
for size in 3
do
for j in 1.0000
do


cat <<END_OF_PROGRAM > L${size}j${j}t${t}h${h}.pbs
#!/bin/sh
# The following line specifies the maximum cpu utilization
#PBS -l cput=1000:00:00
# The following line specifies the maximum cpu utilization
#PBS -l walltime=1000:00:00
# The following line merges the stdout and stderr together
#PBS -j oe
#PBS -l mem=2gb
# The following line specifies the name of the job
#PBS -N ${size}t${t}h${h}
    job=L${size}j${j}t${t}h${h}
cd \${PBS_O_WORKDIR}
#########################################################################

echo \`hostname\` \`date\` >> ./trial\$job
source /opt/intel/cce/9.1.047/bin/iccvars.sh
(time ./\$job.out) >> ./trial\$job
echo \`date\` >> ./trial\$job

cd \${PBS_O_WORKDIR}
if test -e stopfile\$job
then
   date >> ./trial\$job
   echo Chain ended at \`date\` >> ./trial\$job
   rm ./\$job.out
else
/opt/pbs/default/bin/qsub \$job.pbs
fi
exit 0

END_OF_PROGRAM


done
done
done
done
