#! bin/bash
from mpi4py import MPI
import time
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if os.cpu_count() == 48:
   print("\nHello! \n I am your allocated CPU. I belong to the LOGIN NODE. Please, DO NOT USE me to make your computations on int-nano. Use COMPUTE NODES!")
   print("EXAMPLE.  \nrun -p short -n 4 python3 0_2.py \n to allocate and run this script on 4 computing CPUs. \nBye!\n")


if os.cpu_count() != 48:
	if size == 1:
	   print("\nHello! \n I your allocated CPU. I belong to the Compute Node. Please, use us to make your computations on int-nano.")
	   print("EXAMPLE. write: \nrun -p short -n 4 python3 test1.py \n to allocate 4 computing CPUs. \nBye!")
	if size > 1:
	   if rank == 0:
	      print("\nHello! \nI am the master CPU. I do compute nothing, but I have {} slaves (workers)".format(size-1))
	      print("They will now compute 1x1, 2x2, ..., {}x{} for me\n".format((size-1), (size-1)))

	   if rank != 0:
	      time.sleep(1)
	      a = (rank)*(rank)
	      print("I am a worker CPU number {}. I compute {}x{}. \n Result: {}".format(rank, rank, rank, a))


