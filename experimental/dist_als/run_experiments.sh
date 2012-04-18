#!/bin/bash

dist_als=$PWD/../dist_als
dist_fscope_als=$PWD/../dist_fscope_als
dist_sync_als=$PWD/../dist_sync_als



make_hostfile="cat \$PBS_NODEFILE | uniq > machines\n"
hadoop_env="CLASSPATH=\$CLASSPATH:\`hadoop classpath\`\n"
clear_jobs="mpiexec --hostfile machines -n 8 killall dist_als\n 
           mpiexec --hostfile machines -n 8 killall dist_fscope_als\n
           mpiexec --hostfile machines -n 8 killall dist_sync_als\n"

echo -e $header

for d in 10; do
    mkdir dim_$d;
    cd dim_$d;
    for p in $(seq 1 8); do
        mkdir p_$p
        cd p_$p

        mkdir fscope
        cd fscope
        pbs_header="#PBS -d $PWD -N als_fscope_d_${d}_p_${p} -l nodes=${p}:ppn=8 
                    -o bros.ml.cmu.edu:$PWD/output.txt 
                    -e bros.ml.cmu.edu:$PWD/error.txt"
        header="$pbs_header \n $make_hostfile \n $hadoop_env \n $clear_jobs"
        commandstr="mpiexec --hostfile ./machines -n $p
            -x CLASSPATH 
            $dist_fscope_als 
            --matrix hdfs://bros.ml.cmu.edu/data/wikipedia/termdoc_graph 
            --ncpus=8 --tol=0 --maxupdates=2 --dgraphopts=\"(ingress=oblivious)\" 
            --engine=\"(max_clean_fraction=0.5)\" --D $d "
        echo -e $header "\n\n" $commandstr | tee qsub.job
        qsub qsub.job
        cd ..



        mkdir sync
        cd sync
        pbs_header="#PBS -d $PWD -N als_sync_d_${d}_p_${p} -l nodes=${p}:ppn=8 
                    -o bros.ml.cmu.edu:$PWD/output.txt 
                    -e bros.ml.cmu.edu:$PWD/error.txt"
        header="$pbs_header \n $make_hostfile \n $hadoop_env \n $clear_jobs"
        commandstr="mpiexec --hostfile ./machines -n $p
            -x CLASSPATH 
            $dist_sync_als 
            --matrix hdfs://bros.ml.cmu.edu/data/wikipedia/termdoc_graph 
            --ncpus=8 --tol=0 --maxupdates=2 --dgraphopts=\"(ingress=oblivious)\" 
            --engine=\"(max_clean_fraction=0.5)\" --D $d "
        echo -e $header "\n\n" $commandstr | tee qsub.job
        qsub qsub.job
        cd ..


        mkdir locking
        cd locking
        pbs_header="#PBS -d $PWD -N als_locking_d_${d}_p_${p} -l nodes=${p}:ppn=8 
                    -o bros.ml.cmu.edu:$PWD/output.txt 
                    -e bros.ml.cmu.edu:$PWD/error.txt"
        header="$pbs_header \n $make_hostfile \n $hadoop_env \n $clear_jobs"
        commandstr="mpiexec --hostfile ./machines -n $p
            -x CLASSPATH 
            $dist_als 
            --matrix hdfs://bros.ml.cmu.edu/data/wikipedia/termdoc_graph 
            --ncpus=8 --tol=0 --maxupdates=2 --dgraphopts=\"(ingress=oblivious)\" 
            --engine=\"(max_clean_fraction=0.5)\" --D $d "
        echo -e $header "\n\n" $commandstr | tee qsub.job
        qsub qsub.job
        cd ..

        cd ..
    done
    cd ..
done

