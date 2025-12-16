#!/bin/bash
i=1
for file in args*.yaml;
do
echo "$i > $file"
echo
cp submit.sh submit$i.sh
sed -i.bak "s;args.yaml;$file;g" submit$i.sh
sbatch submit$i.sh
((i++))
done
