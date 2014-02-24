#!/bin/bash

# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This script is an example benchmarking of GraphLab for EC2
# for testing scaling
# (C) GraphLab Inc. 2013
# Please send any questions or bug reports to graphlabapi@groups.google.com
# Written by Danny Bickson


############################################################################
# CONFIGURATION
############################################################################
MAX_SLAVES=3  # configure the maximum number of slaves
MAX_RETRY=3   # configure the number of experiemnt repeats
PAGERANK=1    # if 1, runs pagerank
SVD=1         # if 1, runs svd
ALS=1         # if 1, runs als

#It is recommended to define the below two variables for easier setup
#uncomment the below two lines once you set them up
#export AWS_ACCESS_KEY_ID=[ Your access key ]
#export AWS_SECRET_ACCESS_KEY=[ Your access key secret ]
######################################################################

# clean old running instances, if any
echo "y" | ./gl-ec2 -i ~/.ssh/amazonec2.pem -k amazonec2  destroy hpctest  
# launch ec2 cc2.8xlarge image
./gl-ec2 -i ~/.ssh/amazonec2.pem -k amazonec2 -a hpc -s $MAX_SLAVES -t cc2.8xlarge launch hpctest  
# update the GraphLab version to be the latest, recompile, and update slaves
./gl-ec2 -i ~/.ssh/amazonec2.pem -k amazonec2 update hpctest 

# run pagerank benchmarks
if [ $PAGERANK -eq 1 ]; then
for i in `seq 0 1 $MAX_SLAVES`
do
  echo "Running Pagerank"
  for j in `seq 0 1 $MAX_RETRY`
  do
        ./gl-ec2 -i ~/.ssh/amazonec2.pem -k amazonec2 -s $i pagerank_demo hpctest  
  done
done
fi

# run SVD benchmarks
if [ $SVD -eq 1 ]; then
for i in `seq 0 1 $MAX_SLAVES`
do
  echo "Running SVD"
  for j in `seq 0 1 $MAX_RETRY`
  do
        ./gl-ec2 -i ~/.ssh/amazonec2.pem -k amazonec2  -s $i svd_demo hpctest  
  done
done
fi

# run ALS benchmarks
if [ $ALS -eq 1 ]; then
for i in `seq 0 1 $MAX_SLAVES`
do
  echo "Running ALS"
  for j in `seq 0 1 $MAX_RETRY`
  do
     if [ $first_time -eq 1 ]; then
        ./gl-ec2 -i ~/.ssh/amazonec2.pem -k amazonec2  -s $i  als_demo hpctest  
     fi
  done
done
fi

# clean everything
echo "y" | ./gl-ec2 -i ~/.ssh/amazonec2.pem -k amazonec2  destroy hpctest  
