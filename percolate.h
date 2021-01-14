/*
 *  Main header file for percolation code.
 */
#include <mpi.h>

/*
 *  System size L. You can change the value of L to change the problem size.
 */
#define L 432


#define TRUE  1
#define FALSE 0
#define ndims 2


struct neighbors
{
    int up_nbr, down_nbr, lft_nbr, rt_nbr;
};


void cart_topology(int *rank, int size, int *dims, int *periods, struct neighbors *nbrs, MPI_Comm comm, MPI_Comm *comm2D);

void distributing_map_to_old(int M, int N, int old[M+2][N+2], int map[L][L], int rank, int size, int *dims, MPI_Datatype master_vector, MPI_Datatype worker_vector, MPI_Comm comm2D);

void zero_halos(int M, int N, int old[M+2][N+2]);

void swap_halos(int M, int N, int old[M+2][N+2], struct neighbors *nbrs, MPI_Comm comm2D, MPI_Datatype up_down_vector, MPI_Datatype lft_rt_vector);

void collecting_map_from_old(int M, int N, int old[M+2][N+2], int map[L][L], int rank, int size, int *dims, MPI_Datatype master_vector, MPI_Datatype worker_vector, MPI_Comm comm2D);

void test_percolation(int map[L][L]);

/*
 *  Visualisation
 */

void percwrite(char *percfile, int map[L][L], int ncluster);
void percwritedynamic(char *percfile, int **map, int l, int ncluster);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);
