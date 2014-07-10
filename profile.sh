make

if [ $? -ne 0 ]
then
	exit 1
fi


#events="gld_32b,gld_64b,gld_128b,gld_incoherent,gld_coherent,branch,warp_serialize"
#events="gld_32b,gld_64b,gld_128b"

#nvprof --events $events ./smatcher ac -m 8 -p_size 8000 -n 116234496 -alphabet 4

events="gld_incoherent,gld_coherent,branch,warp_serialize"

nvprof --events $events ./smatcher sog -m 8 -p_size 1000 -n 116234496 -alphabet 4


echo ""
echo "gld_32b:  Number of 32 byte global memory load transactions. This increments by 1 for each 32 byte transaction."
echo ""
echo "gld_64b:  Number of 64 byte global memory load transactions. This increments by 1 for each 64 byte transaction."
echo ""
echo "gld_128b:  Number of 128 byte global memory load transactions. This increments by 1 for each 128 byte transaction."
echo ""
echo "gld_incoherent:  Number of non-coalesced global memory loads."
echo ""
echo "gld_coherent:  Number of coalesced global memory loads."
echo ""
echo "branch:  Number of branches taken by threads executing a kernel. This counter will be incremented by one if at least one thread in a warp takes the branch."
echo ""
echo "warp_serialize:  If two addresses of a memory request fall in the same memory bank, there is a bank conflict and the access has to be serialized. This counter gives the number of thread warps that serialize on address conflicts to either shared or constant memory."

