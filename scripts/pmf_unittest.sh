#/bin/sh
cd /tmp/graphlabapi/release/demoapps/pmf
ln -s /mnt/bigbrofs/usr5/graphlab/testdata/netflix .
ln -s /mnt/bigbrofs/usr5/graphlab/testdata/netflixe .
echo "TESTING ALTENATING LEAST SQUARES"
./pmf netflix 0 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false --lambda=0.065
./pmf netflix 1 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false --bptf_burn_in=7
echo "TESTING BPTF"
./pmf netflix 2 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false --bptf_burn_in=7
./pmf netflix 3 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false --bptf_burn_in=7
echo "ALS TENSOR MULT"
./pmf netflix 4 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false --lambda=0.065
echo "TESTING SVD++"
./pmf netflix 5 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false --bptf_burn_in=7
echo "TESTING "
./pmf netflix 6 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false 
./pmf netflix 7 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false 
./pmf netflix 8 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false 
echo "Testing weighted als"
./pmf netflix 9 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false 
echo "Testing sparse factor matrices"
./pmf netflix 10 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false 
./pmf netflix 11 --ncpus=8 --scheduler='round_robin(max_iterations=15,block_size=1)' --float=false 
echo "Testing implicit ratings"
./pmf --unittest=92


