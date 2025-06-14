#!/bin/bash

# Set paths in yaml files for saga
for file in $RGA_LAMBDA_ANALYSIS_HOME/jobs/*/*/*.yaml; do
    sed -i.bak "s;/RGA_LAMBDA_ANALYSIS_HOME;$RGA_LAMBDA_ANALYSIS_HOME;g" $file
    sed -i.bak "s;/RGA_LAMBDA_ANALYSIS_VOL_DIR;$RGA_LAMBDA_ANALYSIS_VOL_DIR;g" $file
done

# Fix for now so that job output directly can be set by saga python libraries
for file in $RGA_LAMBDA_ANALYSIS_HOME/jobs/saga/*/submit.sh; do
    sed -i.bak "s;\$RGA_LAMBDA_ANALYSIS_HOME;$RGA_LAMBDA_ANALYSIS_HOME;g" $file
done
