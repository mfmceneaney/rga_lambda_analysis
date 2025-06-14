#!/bin/bash
cd $RGA_LAMBDA_ANALYSIS_HOME/jobs/c12analysis/mc_rga
i=1
for file in $RGA_MC_DIR/*.hipo;
do
echo "$i > $file"
echo
cp job.sh job$i.sh
cp submit.sh submit$i.sh
sed -i "s;/path/test.hipo;$file;g" job$i.sh
sed -i "s;job.sh;job$i.sh;g" submit$i.sh
sbatch submit$i.sh
((i++))
done
