#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "percolate.h"


/*
 *  This function is to test if percolation occurred by looking for positive numbers
 *  that appear on both the left and right edges
 *
 */

void test_percolation(int map[L][L])
{
    int perc,ilft,irt;

    perc = 0;

    for (ilft=0; ilft < L; ilft++)
    {
      if (map[0][ilft] > 0)
        {
          for (irt=0; irt < L; irt++)
        {
          if (map[0][ilft] == map[L-1][irt])
            {
              perc = 1;
            }
        }
        }
    }

      if (perc != 0)
    {
      printf("percolate: cluster DOES percolate\n");
    }
      else
    {
      printf("percolate: cluster DOES NOT percolate\n");
    }
}