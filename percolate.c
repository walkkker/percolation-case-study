/*
 * 2D decomposition parallel program to test for percolation of a cluster
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "percolate.h"



int main(int argc, char *argv[])
{

  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */

  int i, j, nhole, step, maxstep, oldval, newval;
  int nchangelocal, nchange, printfreq;
  long sum, local_sum; // Used to calculate the average value of the map array
  double r;


 /*
  * variables for timing program to measure the parallel performance
  */
  double tstart, tstop;

  /*
   *  MPI variables
   */

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm comm2D;

  int size, rank;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);


/**************************Cartesian topology*********************************/
  
  int periods[ndims], dims[ndims];
  struct neighbors *nbrs, nbrs1; nbrs = &nbrs1;

  /*
   *  Put all processes into a Cartesian topology and get the number of processes in each dimension.
   *  This function also computes the 4 neighbours for each process, where MPI_PROC_NULL is assigned automatically.
   */

  cart_topology(&rank, size, dims, periods, nbrs, comm, &comm2D);

/*****************************************************************************/


   /*
    * Check if L is an exact multiple of the number of processes  
    * to ensure our program can run correctly.
    */

    if (L%dims[0]!=0 || L%dims[1]!=0)
    {
      if (rank == 0)
      {
        printf("percolate: ERROR L is not an exact multiple of the number of processes\n");
      }

      MPI_Finalize();
      return 0;
    }


    /*
     *  If the number of processes, P, can evenly divide the grid size, L. 
     *  Calculate the length and width of divided rectangulars.
     */

    int M = L/dims[0];
    int N = L/dims[1];

    /*
    *  Define the main arrays for the simulation
    */

    int old[M+2][N+2], new[M+2][N+2];

    int map[L][L];



    if (argc != 2)
    {
      if (rank == 0)
    {
      printf("Usage: percolate <seed>\n");
    }

      MPI_Finalize();
      return 0;
    }

  /*
   *  Update for a fixed number of iterations
   */

  maxstep = 4*L;
  printfreq = 100;

  if (rank == 0)
    {
      printf("percolate: running on %d process(es)\n", size);

      /*
       *  Set most important value: the rock density rho (between 0 and 1)
       */

      rho = 0.411;

      /*
       *  Set the randum number seed and initialise the generator
       */

      seed = atoi(argv[1]);

      printf("percolate: L = %d, rho = %f, seed = %d, maxstep = %d\n",
         L, rho, seed, maxstep);

      rinit(seed);

     /*
      *  Initialise map with density rho. Zero indicates rock, a positive
      *  value indicates a hole. For the algorithm to work, all the holes
      *  must be initialised with a unique integer
      */

      nhole = 0;

      for (i=0; i < L; i++)
        {
          for (j=0; j < L; j++)
            {
              r=uni();
          
              if(r < rho)
            {
              map[i][j] = 0;
            }
              else
            {
              nhole++;
              map[i][j] = nhole;
            }
            }
        }

      printf("percolate: rho = %f, actual density = %f\n",
          rho, 1.0 - ((double) nhole)/((double) L*L) );
    }



   /*
    *  Master_vector is used for map array, while worker_vector is used for old array.
    *  Through these two vectors, we can transfer data directly between map and old.
    */

    MPI_Datatype master_vector, worker_vector;
    int blocklength = N;
    int stride = L;
    int count = M;
    MPI_Type_vector (count, blocklength, stride, MPI_INT, &master_vector);
    MPI_Type_commit(&master_vector);

    blocklength = N;
    stride = N+2;
    count = M;
    MPI_Type_vector (count, blocklength, stride, MPI_INT, &worker_vector);
    MPI_Type_commit(&worker_vector);
    



/**********Distributing map:transfer data directly from map to the internal region of old***********/

    distributing_map_to_old(M, N, old, map, rank, size, dims, master_vector, worker_vector, comm2D);

/************************************Distribution finished******************************************/

    step = 1;
    nchange = 1;
    sum =1;

    /* up-down vector. Used for changing data with up and down neighbors */
    MPI_Datatype up_down_vector;
    blocklength = 1;
    stride = N+2;
    count = M;
    MPI_Type_vector (count, blocklength, stride, MPI_INT, &up_down_vector);
    MPI_Type_commit(&up_down_vector);

    /* left-right vector. Used for changing data with left and right neighbors */
    MPI_Datatype lft_rt_vector;
    blocklength = N;
    stride = N;
    count = 1;
    MPI_Type_vector (count, blocklength, stride, MPI_INT, &lft_rt_vector);
    MPI_Type_commit(&lft_rt_vector);


   /*
    * Line up at the start line and fire the gun and start the clock.
    */
    MPI_Barrier(comm2D);
    tstart = MPI_Wtime();


  while (step <= maxstep)
    {

/*************************************Swap halos*********************************************/

      swap_halos(M, N, old, nbrs, comm2D, up_down_vector, lft_rt_vector);

/**********************************Halo swaps finished***************************************/

    local_sum = 0;

    nchangelocal = 0;

    for (i=1; i<=M; i++)
    {
      for (j=1; j<=N; j++)
        {
          oldval = old[i][j];
          newval = oldval;

        /*
         * Set new[i][j] to be the maximum value of old[i][j]
         * and its four nearest neighbours
         */

          if (oldval != 0)
        {
          if (old[i-1][j] > newval) newval = old[i-1][j];
          if (old[i+1][j] > newval) newval = old[i+1][j];
          if (old[i][j-1] > newval) newval = old[i][j-1];
          if (old[i][j+1] > newval) newval = old[i][j+1];

          if (newval != oldval)
            {
              ++nchangelocal;
            }
        }

          new[i][j] = newval;

          local_sum += new[i][j];
        }
    }

    /*
     *  Compute global number of changes to all processes
     */

     MPI_Allreduce(&nchangelocal, &nchange, 1, MPI_INT, MPI_SUM, comm2D);

    /*
     *  Calculate the sum of the map array. average = sum/(L*L)
     */

     MPI_Reduce(&local_sum, &sum, 1, MPI_LONG, MPI_SUM, 0, comm2D);

    /*
     *  Report progress every now and then
     */

    if (step % printfreq == 0)
    {
      if (rank == 0)
        {
          printf("percolate: number of changes on step %d is %d\n",step, nchange);
          printf("percolate: the average value of the map array on step %d is %.2lf\n", step, ((double)sum/(L*L)));
        }
    }

     /*
      *  Copy back in preparation for next step, omitting halos
      */

    for (i=1; i<=M; i++)
    {
      for (j=1; j<=N; j++)
        {
          old[i][j] = new[i][j];
        }
    }

   /*
    * Set stopping criterion.
    * Computation stops as soon as there are no changes, rather than doing a fixed number of steps.
    */

    if(nchange == 0)
    {
      break;
    }
    
    step++; 

    }



   /*
    * Wait for every process to finish and stop the clock.
    */
    MPI_Barrier(comm2D);
    tstop = MPI_Wtime();


  /*
   *  We set a maximum number of steps to ensure the algorithm always
   *  terminates. However, if we hit this limit before the algorithm
   *  has finished then there must have been a problem (e.g. maxstep
   *  is too small)
   */

  if (rank == 0)
    {
      if (nchange != 0)
    {
      printf("percolate: WARNING, steps = %d reached before nchange = 0\n",
         maxstep);
    }
    }



/*********Collecting map: transfer data directly from the internal region of old to map**********/

    collecting_map_from_old(M, N, old, map, rank, size, dims, master_vector, worker_vector, comm2D);

/*****************************************Collection finished************************************/




  if (rank == 0)
    {

     /*
      *  Test to see if percolation occurred by looking for positive numbers
      *  that appear on both the left and right edges
      */

      test_percolation(map);


     /*
      *  Write the map to the file "map.pgm", displaying the two
      *  largest clusters. If the last argument here was 3, it would
      *  display the three largest clusters etc. The picture looks
      *  cleanest with only a single cluster, but multiple clusters
      *  are useful for debugging.
      */

      percwrite("map.pgm", map, 3);

      printf("percolate: Execution time between distribution and collection (not included) is %f\n", (float)(tstop-tstart));
    }

  MPI_Finalize();

  return 0;
}
