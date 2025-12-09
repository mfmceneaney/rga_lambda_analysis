#!/bin/bash

# Got to job directory
cd $RGA_LAMBDA_ANALYSIS_HOME/jobs/simulation/mc_rga_sss

# Loop root files
j=0
for file in $RGA_LAMBDA_ANALYSIS_ROOT_FILES_DIR/*.root;
do
offset=$(($RGA_LAMBDA_ANALYSIS_NFILES * $j))
# Loop number of lund files
i=1
while [ $i -le $RGA_LAMBDA_ANALYSIS_NFILES ];
do
idx=$(($offset + $i))
echo "$idx > $PWD/submit$idx.sh"
echo

# Copy scripts
cp job.sh job$idx.sh
cp submit.sh submit$idx.sh

# Replace variables
sed -i "s;INFILE=\"file.root\";INFILE=$file;g" job$idx.sh
sed -i "s;MCINDEX=0;MCINDEX=$i;g" job$idx.sh
sed -i "s;PREFIX=\"out_\";PREFIX=\"out_${j}_\";g" job$idx.sh
sed -i "s;job.sh;job$idx.sh;g" submit$idx.sh

# Submit job
sbatch submit$idx.sh
((i++))
done
((j++))
done
