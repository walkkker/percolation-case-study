#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "percolate.h"


/*
 *  This function transfers data directly from the internal regions of old of all processes to the 'map' array of controller (rank 0),
 *  so it can achieve the gather function.
 */

void collecting_map_from_old(int M, int N, int old[M+2][N+2], int map[L][L], int rank, int size, int *dims, MPI_Datatype master_vector, MPI_Datatype worker_vector, MPI_Comm comm2D)
{
    int i,j;
    int coords[2];
    int cart_rank;


    MPI_Request request;
    MPI_Request array_of_requests[size];
    MPI_Status status;
    MPI_Status array_of_statuses[size];
    
    /* All processes send the internal region of old to the controller (rank 0) */
    MPI_Issend ( &old[1][1], 1, worker_vector,0, 0, comm2D,&request);
    
    
   /*
    *  The controller receives the internal regions of 'old' arrays from different processes, and 
    *  put them in their corresponding positions in 'map' according to the Cartesian topology.
    */
    if (rank == 0)
    {
        for(i=0; i<dims[0]; i++)
        {
            for(j=0; j<dims[1]; j++)
            {
                coords[0] = i;
                coords[1] = j;
                MPI_Cart_rank(comm2D, coords, &cart_rank);
                MPI_Irecv (&map[i*M][j*N], 1, master_vector,cart_rank, 0, comm2D,&array_of_requests[cart_rank]);
            }
        }
    }

    MPI_Wait(&request, &status);
    
    if(rank == 0){
        MPI_Waitall (size, array_of_requests,array_of_statuses);
    }
}