#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "percolate.h"

/*
 *  Use non-blocking communications with each process's four neighbors in Cartesian topology to implement halo swaps.
 */

void swap_halos(int M, int N, int old[M+2][N+2], struct neighbors *nbrs, MPI_Comm comm2D, MPI_Datatype up_down_vector, MPI_Datatype lft_rt_vector)
{

    MPI_Request swap_requests[8];
    MPI_Status swap_statuses[8];

    /* Non-blocking communications with left neighbor */
    MPI_Issend(&old[1][1],1,lft_rt_vector,nbrs->lft_nbr,0,comm2D,&swap_requests[0]); 
    MPI_Irecv(&old[0][1],1,lft_rt_vector,nbrs->lft_nbr,1,comm2D,&swap_requests[1]);

    /* Non-blocking communications with right neighbor */
    MPI_Issend(&old[M][1],1,lft_rt_vector,nbrs->rt_nbr,1,comm2D,&swap_requests[2]); 
    MPI_Irecv(&old[M+1][1],1,lft_rt_vector,nbrs->rt_nbr,0,comm2D,&swap_requests[3]);

    /* Non-blocking communications with up neighbor */
    MPI_Issend(&old[1][N],1,up_down_vector,nbrs->up_nbr,2,comm2D,&swap_requests[4]);
    MPI_Irecv(&old[1][N+1],1,up_down_vector,nbrs->up_nbr,3,comm2D,&swap_requests[5]);

    /* Non-blocking communications with down neighbor */
    MPI_Issend(&old[1][1],1,up_down_vector,nbrs->down_nbr,3,comm2D,&swap_requests[6]);
    MPI_Irecv(&old[1][0],1,up_down_vector,nbrs->down_nbr,2,comm2D,&swap_requests[7]);

    MPI_Waitall(8, swap_requests, swap_statuses);
}