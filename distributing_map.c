#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "percolate.h"


/*
 * The controller distributing the divided map to all processes. Data is transfered from map to old directly.
 */

void distributing_map_to_old(int M, int N, int old[M+2][N+2], int map[L][L], int rank, int size, int *dims, MPI_Datatype master_vector, MPI_Datatype worker_vector, MPI_Comm comm2D)
{
    MPI_Request request;
    MPI_Status status;

    MPI_Request array_of_requests[size];
    MPI_Status array_of_statuses[size];

    int coords[2];
    int cart_rank;
	int i,j;

    if(rank == 0){
        for(i=0; i<dims[0]; i++){
            for(j=0; j<dims[1]; j++){
                coords[0]=i;
                coords[1]=j;
                MPI_Cart_rank(comm2D, coords, &cart_rank);
                MPI_Issend ( &map[i*M][j*N], 1, master_vector,cart_rank, 0, comm2D,&array_of_requests[cart_rank]);
            }
        }
    }
    MPI_Irecv (&old[1][1], 1, worker_vector,0, 0, comm2D,&request);


  /*
   * Overlapping communication and calculation. Doing Communication while set the halo values to zero.
   */
    zero_halos(M, N, old);
    

    if(rank == 0)
    {
        MPI_Waitall (size, array_of_requests, array_of_statuses);
    }
    MPI_Wait(&request, &status);
}


/*
 *  Zero the halos of "old" array.
 */

void zero_halos(int M, int N, int old[M+2][N+2])
{
   int i,j;
   
   for (i=0; i <= M+1; i++)  // zero the bottom and top halos
    {
      old[i][0]   = 0;
      old[i][N+1] = 0;
    }

   for (j=0; j <= N+1; j++)  // zero the left and right halos
    {
      old[0][j]   = 0;
      old[M+1][j] = 0;
    }
}