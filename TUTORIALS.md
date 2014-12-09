# GraphLab PowerGraph Tutorials

##Table of Contents
* [Deploying on AWS EC2 Cluster](#ec2)
* [Deploying in a Cluster](#cluster)
* [Deploying on a single multicore machine](#multicore)
* [Benchmarking on AWS EC2](#benchmarking)
* [Fine tuning GraphLab PowerGraph performance](#perf_tuning)

<a name="ec2"></a>
# Deploying in AWS EC2 Cluster

## Step 0: Requirements
* You should have Amazon EC2 account eligible to run on us-east-1a zone.

* Find out using the Amazon AWS console your AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY (under your account name on the top right corner-&gt; security credentials -&gt; access keys)

* You should have a keypair attached to the zone you are running on (in our example us-east-1a) as explained <a {{ trackClick() }} href="http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html">here</a>. You will need to know your keypair name (graphlabkey in our example), and the location of the private key (~/.ssh/graphlabkey.pem in our example).

* Install [boto](https://pypi.python.org/pypi/boto/). This is the AWS Python client. To install, run: 

```
sudo pip install boto
```

* Download and install GraphLab PowerGraph using the instructions in the [README.md](README.md).


## Step 1: Environment Setup

Edit your .bashrc or .bash_profile or .profile files (remember to source it after editing, using the bash command “source &lt;filename&gt;”)

```
export AWS_ACCESS_KEY_ID=[ Your access key ]
export AWS_SECRET_ACCESS_KEY=[ Your access key secret ]
```

## Step 2: Start the cluster

```
cd ~/graphlabapi/scripts/ec2
./gl-ec2 -i ~/.ssh/graphlab.pem -k graphlabkey  -s 1 launch launchtest
```

(In the above command, we created a 2-node cluster in us-east-1a zone. -s is the number of slaves, launch is the action, and launchtest is the name of the cluster)

## Step 3: Update GraphLab PowerGraph

```
./gl-ec2 -i ~/.ssh/graphlab.pem -k graphlabkey update launchtest
```

## Step 4: Run Alternating Least Squares Demo

This step runs ALS (alternating least squares) in a cluster using small netflix subset.
It first downloads the data from the web: [http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.train](http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.train) and [http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.validate](http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.validate), copy it into HDFS, and runs 5 alternating least squares iterations:

```
./gl-ec2 -i ~/.ssh/graphlab.pem -k graphlabkey als_demo launchtest
```

After the run is completed, login to the master node and view the output files in the folder ~/graphlabapi/release/toolkits/collaborative_filtering/ The algorithm and exact format is explained in the API docs.

## Step 5: Shutdown the Cluster

```
./gl-ec2 -i ~/.ssh/graphlab.pem -k grpahlabkey destroy launchtest
```

## Other Useful Commands:

Login into the master node using

```
./gl-ec2 -i ~/.ssh/graphlab.pem -s 1 login launchtest
```

<a name="cluster"></a>
# Deploying in a Cluster

## Step 0: Install GraphLab PowerGraph on one of your cluster nodes.

Install GraphLab PowerGraph, using instructions in the [README.md](README.md), on your master node (one of your cluster machines).

## Step 1: Copy GraphLab PowerGraph files to all machines.

1) Create a file called in your home directory called “machines” with the names of all the MPI nodes participate in the computation.

For example:

```
cat ~/machines
mynode1.some.random.domain
mynode2.some.random.domain
...
mynode18.some.random.domain
```
2) Verify you have the machines files from section 1) in your root folder of all of the machines.

3) You will need to setup password-less SSH between the master node and all other machines.

Verify it is possible to ssh without password between any pairs of machines. These [instructions](http://www.linuxproblem.org/art_9.html) explain how to setup ssh without passswords.

Before proceeding, verify that this is setup correctly; check that the following connects to the remote machine without prompting for a password:

```
# from machine mynode1.some.random.domain
ssh mynode2.some.random.domain
```

4) On the node you installed GraphLab on, run the following commands to copy GraphLab files to the rest of the machines:

```
cd ~/graphlab/release/toolkits
~/graphlab/scripts/mpirsync
cd ~/graphlab/deps/local
~/graphlab/scripts/mpirsync
```

This step will only work if the file you created in step 1 was named "machines" and located in your home directory.

In order for mpirsync to run properly all machines must have all network ports open.

## Step 2a: Run PageRank on a synthetic graph

This step runs the [PageRank](http://en.wikipedia.org/wiki/PageRank) algorithm on a synthetic generated graph of 100,000 nodes. It spawns two GraphLab mpi instances (-n 2).
```
mpiexec -n 2 -hostfile ~/machines /path/to/pagerank --powerlaw=100000
```

## Step 2: Run GraphLab PowerGraph ALS using subset of Netflix data

This step runs ALS (alternating least squares) in a cluster using small netflix susbset.
It first downloads an anonymized, synthetic Netflix dataset from the web: [http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.train](http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.train) and [http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.validate](http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.validate), and runs 5 alternating least squares iterations. After the run is completed, you can login into any of the nodes and view the output files in the folder ~/graphlab/release/toolkits/collaborative_filtering/

 ```
 cd /some/ns/folder/
mkdir smallnetflix
cd smallnetflix/
wget http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.train
wget http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.validate
```
Now run GraphLab:

````
mpiexec -n 2 -hostfile ~/machines /path/to/als  --matrix /some/ns/folder/smallnetflix/ --max_iter=3 --ncpus=1 --minval=1 --maxval=5 --predictions=out_file
```
Where -n is the number of MPI nodes, and –ncpus is the number of deployed cores on each MPI node.

machines is a file which includes a list of the machines you like to deploy on (each machine in a new line)

Note: this section assumes you have a network storage (ns) folder where the input can be stored.
Alternatively, you can split the input into several disjoint files, and store the subsets on the cluster machines.

Note: Don’t forget to change /path/to/als and /some/ns/folder to your actual folder path!

Note: For mpich2, use -f instead of -hostfile.

## Step 3:

[Fine tuning graphlab deployment](#perf_tuning).

## Errors and their resolution:

### Error:

```
/mnt/info/home/daroczyb/als: error while loading shared libraries: libevent_pthreads-2.0.so.5: cannot open shared object file: No such file or directory
```

**Solution:**

You should define LD_LIBRARY_PATH to point to the location of libevent_pthreads, this is done with the -x mpi command, for example:

```
mpiexec --hostfile machines -x LD_LIBRARY_PATH=/home/daroczyb/graphlab/deps/local/lib/ /mnt/info/home/daroczyb/als /mnt/info/home/daroczyb/smallnetflix_mm.train
```

### Error:

```
mnt/info/home/daroczyb/als: error while loading shared libraries: libjvm.so: cannot open shared object file: No such file or directory
```

**Solution:**

Point LD_LIBRARY_PATH to the location of libjvm.so using the -x mpi command:

```
mpiexec --hostfile machines -x LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/daroczyb/graphlab/deps/local/lib/:/usr/lib/jvm/java-7-openjdk-amd64/jre/lib/amd64/server/ /mnt/info/home/daroczyb/als /mnt/info/home/daroczyb/smallnetflix_mm.train
```

### Error:

```
problem with execution of /graphlab/release/toolkits/collaborative_filtering/als  on  debian1:  [Errno 2] No such file or directory
```

**Solution:**

You should verify the executable is found on the same path on all machines.

### Error:

a prompt asking for password when running mpiexec

**Solution:** Use the following [instructions](http://www.linuxproblem.org/art_9.html) to allow connection with a public/private key pair (no password).

### Error:

```
Exception in thread "main" java.lang.IllegalArgumentException: Wrong FS: hdfs://[domain]:9000/user/[user_name]/data.txt, expected: file:///
    at org.apache.hadoop.fs.FileSystem.checkPath(FileSystem.java:381)
    at org.apache.hadoop.fs.RawLocalFileSystem.pathToFile(RawLocalFileSystem.java:55)
    at org.apache.hadoop.fs.RawLocalFileSystem.listStatus(RawLocalFileSystem.java:307)
    at org.apache.hadoop.fs.FileSystem.listStatus(FileSystem.java:842)
    at org.apache.hadoop.fs.FileSystem.listStatus(FileSystem.java:867)
    at org.apache.hadoop.fs.ChecksumFileSystem.listStatus(ChecksumFileSystem.java:487)
    Call to org.apache.hadoop.fs.FileSystem::listStatus failed!
    WARNING: distributed_graph.hpp(load_from_hdfs:1889): No files found matching hdfs://[domain]:9000/user/[user_name]/data.txt
```

**Solution:**
Verify you classpath includes all hadoop required folders.

### Error:

Just after TCP Communication layer is constructed: 
```
BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES, EXITCODE: 11, CLEANING UP REMAINING PROCESSES, YOU CAN IGNORE THE BELOW CLEANUP MESSAGES
```
or:

```
[xyzserver:22296] *** Process received signal *** mpiexec noticed that process rank 0 with PID 22296 on node xyzserver exited on signal 11 (Segmentation fault).
```

**Solution:**

Check that all machines have access to, or are using the same binary

<a id="multicore"></a>
#Deployment on a single multicore machine

## Preliminaries:

## Step 0: Install GraphLab on one of your cluster nodes.

Using the instructions [here](/projects/source.html) on your master node (one of your cluster machines), except invoke the  configure script with the ‘–no_mpi’ flag.
Don’t forget to use
```
./configure --no_mpi
```

when configuring GraphLab.

## Step 1: Run GraphLab ALS

This step runs ALS (alternating least squares) in a cluster using small netflix susbset. It first downloads the data from the web, runs 5 alternating least squares iterations. After the run is completed, the output files will be created in the running folder (the folder graphlab/release/toolkits/collaborative_filtering/) 

```
cd graphlab/release/toolkits/collaborative_filtering/
mkdir smallnetflix
cd smallnetflix/
wget http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.train
wget http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.validate
cd ..
```

Now run GraphLab:

```
./als --matrix ./smallnetflix/ --max_iter=5 --ncpus=1 --predictions=out_file
```
    
Where –ncpus is the number of deployed cores.

<a id="benchmarking"></a>
# Benchmarking on AWS EC2

A commonly repeating task is evaluation of GraphLab performance and scaling properties on a cluster. To help jump start benchmarking we have created this tutorial.

## Step 0: Requirements

1. You should have Amazon EC2 account eligible to run on us-west zone.
2. Find out using the Amazon AWS console your AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY (under your account name on the top right corner-> security credentials -> access keys)
3. You should have a keypair attached to the zone you are running on (in our example us-west) as explained [here](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html). You will need to know your keypair name (amazonec2 in our example), and the location of the private key (~/.ssh/amazonec2.pem in our example).
4. Install boto. This is the AWS Python client. To install, run: `sudo pip boto`.
5. Download and install GraphLab  using the instructions [here](/projects/source.html).

## Step 1: Recommended setting

We recommend using high performance computing instances (like cc2.8xlarge) since we observed a significant improved performance especially related to variation in cluster load and network utilization. The scripts also allow using regular instances.

To avoid ec2 unexpected loads, we recommend repeating each experiment a few times and computing the average.

## Step 2: Environment Setup

Edit your .bashrc or .bash_profile or .profile files (remember to source it after editing, using the bash command “source <filename>”)

```
export AWS_ACCESS_KEY_ID=[ Your access key ]
export AWS_SECRET_ACCESS_KEY=[ Your access key secret ]
```

## Step 3: configure benchmarking

Edit the [benchmark_ec2.sh](https://github.com/graphlab-code/graphlab/blob/master/scripts/ec2/benchmark_ec2.sh) script found under graphlab/scripts/ec2
1. Select the requested algorithms of the following options:
 ```
ALS=1 # alternating least squares
SVD=1 # singular value decomposition
PR=1  # pagerank
```
(Setting an algorithm to 0 will disable its run).
2. Select the number of slaves (any number between 0 to n) by setting the MAX_SLAVES variable.
3. Select the number of experiment repeats (any number between 0 to n) by setting the MAX_RETRY variable. The benchmarking script, spawns an ec2 cluster of size n machines, and then tests the requested algorithm using 0, 1, … n-1 slaves. Each experiment is repeated MAX_RETRY times.

### Step 3: Perform benchmarking

```
cd ~/graphlabapi/scripts/ec2
./benchmark_ec2.sh
```
It is advised to redirect the benchmarking output to file, for example on bash:

 ```
 ./benchmark_ec2 > output 2>&1
 ```

### Step 4: Processing the results

For detecting final runtime for ALS/SVD

```
grep "Runtime" output
```
For detecting final runtime for PR:

```
grep "Finished Running" output
```
You will need to manually compute the average runtime for each case. A recommended metric to use is the “speedup” curve, which is the time for executing on a single machine divided by the time executing on k machines. The optimal result is linear speedup, namely running on k machines speeds up the algorithm k times vs. running on a single machine.

### Step 5: behind the scenes

Here is a more detailed explanation of the benchmarking process. The benchmarking is calling gl-ec2 script which calls [gl_ec2.py](https://github.com/graphlab-code/graphlab/blob/master/scripts/ec2/gl_ec2.py) script.
1. The “launch” command to start a graphlab cluster with X machines.
2. The “update” command to get the latest version of graphlab from git, recompile it, and disseminate the binary to the salves
3. The “als_demo”, “svd_demo”, “pagerank_demo” command benchmark ALS/SVD/PR algorithms. It first downloads a dataset from the web and then calls graphlab with the right command lines to issue a run on the downloaded dataset. For PR we use the [LiveJournal](http://snap.stanford.edu/data/soc-LiveJournal1.html) dataset. For ALS/SVD we use a [netflix like synthetic sample](http://www.select.cs.cmu.edu/code/graphlab/datasets/smallnetflix_mm.train).
4. In case you would like to benchmark a different dataset, you can edit the dataset URL in the gl_ec2.py example.
5. In case you would like to benchmark a different algorithm, you can add an additional youralgo_demo section into the gl_ec2.py script.
6. In case you would like to bechmark a regular instance, simply change the following line in gl_ec2.py from

````
./gl-ec2 -i ~/.ssh/amazonec2.pem -k amazonec2 -a hpc -s $MAX_SLAVES -t cc2.8xlarge launch hpctest
```
to:
```
./gl-ec2 -i ~/.ssh/amazonec2.pem -k amazonec2  -s $MAX_SLAVES -t m1.xlarge launch hpctest
```

### Advanced topics.

In case you like to work in a different ec2 region (than the default us-west):

For us-east region, those are the provided AMIs:

Standard: ami-31360458, high performance: ami-39360450.

You should

1. add the following line just before: [gl_ec2.py](https://github.com/graphlab-code/graphlab/blob/master/scripts/ec2/gl_ec2.py#L223)
    
```
opts.ami = "ami-31360458"
```
2. run with the additional command line argument:
```
-r us-east-1
```

### Support

If you encounter any problem when trying to run this benchmarking feel free to post on [forum.graphlab.com](http://forum.graphlab.com)

<a id="perf_tuning"></a>
# Fine tuning GraphLab PowerGraph performance

This section contains tips and examples on how to setup GraphLab properly on your cluster and how to squeeze performance.

## 0: Compile in release

Verify you compiled graphlab in the release subfolder (and not in debug subfolder). Compiling in release may speed execution up to x10 times!

Tip: Always compile in release when testing performance.

## 1: Understanding input graph loading

GraphLab PowerGraph has built in parallel loading of the input graph. However, for efficient parallel loading, the input file should be split into multiple disjoint sub files. When using a single input file, the graph loading becomes serial (which is bad!).

Each MPIinstance has a single loader of the input graph attached to it (does not matter how many cpus are used by that MPI instance).

Tip: Always split your input file into at least as many MPI processes you are using.

## 2: Verify MPI is working correctly

You can test your MPI setup as follows:

1. Compile the release/demoapps/rpc subfolder (using “cd release/demoapps/rpc/; make”). Copy the files generated by the compile to all machines.
2. Run:

```
mpiexec -n 2 --hostfile ~/machines  /home/ubuntu/graphlab/release/demoapps/rpc/rpc_example1
```
As part of the output, you should see something like this:

```
TCP Communication layer constructed.
TCP Communication layer constructed.

10
5 plus 1 is : 6
11 plus 1 is : 12
```

If you get something else, please report an error as explained below

## 3: Fine tuning of the partitioning.

Previous to the program execution, the graph is first loaded into memory and partitioned into the different cluster machines. It is possible to try different partitioning strategies. This is done using the following flags:

```
--graph_opts="ingress=oblivious
```

or

````
--graph_opts="ingress=grid" # works for power of 2 sized cluster i.e. 2,4,8,.. machines
```

For different graphs, different partitioning methods may give different performance gains.

## 4: Setting ncpus

The –ncpus option let you set the number of cores used to perform computation. Prior to 2.1.4644 this defaults to 2. After 2.1.4644, this defaults to #cores – 2. When run in the distributed setting, the maximum number this should be set to is #cores – 2 since 2 cores should be reserved for communication.
