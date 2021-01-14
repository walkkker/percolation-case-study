#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "percolate.h"

/*
 *  Create a Cartesian topology for all processes and get the number of processes in each dimension.
 *  This function also computes the 4 neighbours for each process, where MPI_PROC_NULL is assigned automatically.
 */

void cart_topology(int *rank, int size, int *dims, int *periods, struct neighbors *nbrs, MPI_Comm comm, MPI_Comm *comm2D)
{
   /*
    * Cartesian topology
    */

    dims[0] = 0;
    dims[1] = 0;
    periods[0] = FALSE;       // Non-periodic
    periods[1] = TRUE;        // Periodic boundary conditions in the vertical direction.
    int reorder = FALSE;      // FALSE
    int disp = 1;             // Shift by 1

  /*
   * Computing neighbours with a Cartesian topology
   * and MPI_Cart_shift, where MPI_PROC_NULL is assigned automatically.
   */

    MPI_Dims_create(size,ndims,dims);
    MPI_Cart_create(comm,ndims,dims,periods,reorder,comm2D);
    MPI_Comm_rank(*comm2D, rank);
    MPI_Cart_shift(*comm2D,0,disp,&nbrs->lft_nbr,&nbrs->rt_nbr);
    MPI_Cart_shift(*comm2D,1,disp,&nbrs->down_nbr,&nbrs->up_nbr);
}