#include <stdio.h>
#include <stdlib.h>

#include "percolate.h"

/*
 *  Function to write a percolation map in greyscale Portable Grey Map
 *  (PGM) format. The largest "ncluster" clusters are identified and
 *  shown as shades of grey against a black background, with the
 *  largest cluster shown in white.
 */

#define MAXNCLUSTER 9          // Must be able to identify by a single digit
int foundcluster[MAXNCLUSTER];

#define MAXCLUSTERID (L*L)    // Extreme case of zero density
int clustersize[MAXCLUSTERID+1];

void percwrite(char *percfile, int map[L][L], int ncluster)
{
  FILE *fp;

  int i, j, colour, npix;
  int clusterid, icluster, maxcluster, prevcluster;

  static int pixperline = 32; // PGM format limits to 70 characters per line

  if (ncluster > MAXNCLUSTER)
    {
      printf("percwrite: WARNING ncluster too large, resetting to %d\n",
	     MAXNCLUSTER);

      ncluster = MAXNCLUSTER;
    }

  if (ncluster > 1)
    {
      printf("percwrite: visualising the largest %d clusters\n", ncluster);
    }
  else
    {
      printf("percwrite: only visualising the largest cluster\n");
    }

  /*
   * Count up the size of each cluster
   */

  for (i=0; i<=MAXCLUSTERID; i++)
    {
      clustersize[i] = 0;
    }

  for (i=0; i < L; i++)
    {
      for (j=0; j < L; j++)
	{
	  clusterid = map[i][j];

	  if (clusterid > 0)
	    {
	      clustersize[clusterid]++;
	    }
	}
    }

  /*
   * Find the size of the "ncluster" largest clusters (by brute force!)
   */

  prevcluster = MAXCLUSTERID+1; // Larger than the largest possible cluster id

  for (icluster = 0; icluster < ncluster; icluster++)
    {
      maxcluster = 0;

      for (i=0; i<=MAXCLUSTERID; i++)
	{
	  if (clustersize[i] > maxcluster && clustersize[i] < prevcluster)
	    {
	      maxcluster = clustersize[i];
	    }
	}

      foundcluster[icluster] = maxcluster;
      prevcluster = maxcluster;
    }

  if (ncluster > 1)
    {
      printf("percwrite: cluster sizes are ");
    }
  else
    {
      printf("percwrite: maximum cluster size is ");
    }

  for (icluster = 0; icluster < ncluster-1; icluster++)
    {
      printf("%d, ", foundcluster[icluster]);
    }
  printf("%d\n", foundcluster[ncluster-1]);

  /*
   *  Write the file
   */

  printf("percwrite: opening file <%s>\n", percfile);

  fp = fopen(percfile, "w");

  printf("percwrite: writing data ...\n");

  /*
   *  Start with the PGM header
   */

  fprintf(fp, "P2\n");
  fprintf(fp, "%d %d\n%d\n", L, L, ncluster);

  /*
   *  Now write the cells to file so that map[0][0] is in the
   *  bottom-left-hand corner and map[l-1][l-1] is in the
   *  top-right-hand corner
   */

  npix = 0;

  for (j=L-1; j >= 0; j--)
    {
      for (i=0; i < L; i++)
	{
	  clusterid = map[i][j];

	  /*
	   * Write out the largest cluster(s), shading appropriately
	   */

	  colour = 0;

	  if (clusterid > 0)
	    {
	      for (icluster = 0; icluster < ncluster; icluster++)
		{
		  if (clustersize[clusterid] == foundcluster[icluster])
		    {
		      // Largest (first) cluster is white

		      colour = ncluster - icluster;
		    }
		}
	    }

	  npix++;

	  // Make sure lines wrap after "npix" pixels

	  if (npix == 1)
	    {
	      fprintf(fp, "%1d", colour);
	    }
	  else if (npix < pixperline)
	    {
	      fprintf(fp, " %1d", colour);
	    }
	  else
	    {
	      fprintf(fp, " %1d\n", colour);
	      npix = 0;
	    }
	}
    }

  if (npix != 0) fprintf(fp, "\n");

  printf("percwrite: ... done\n");

  fclose(fp);
  printf("percwrite: file closed\n");
}


/*
 *  Function to write a percolation map in greyscale Portable Grey Map
 *  (PGM) format. The largest "ncluster" clusters are identified and
 *  shown as shades of grey against a black background, with the
 *  largest cluster shown in white.
 *
 *  Note that this version expects the map array to have been
 *  dynamically allocated, e.g. using the arralloc() routine:
 *
 *  int **map;
 *  map = (int **) arralloc(sizeof(int), 2, L, L);
 *  ...
 *  percwritedynamic("map.pgm", map, L, 1);
 */

#define MAXNCLUSTER 9          // Must be able to identify by a single digit
int foundcluster[MAXNCLUSTER];

void percwritedynamic(char *percfile, int **map, int l, int ncluster)
{
  FILE *fp;

  int i, j, colour, npix;
  int clusterid, icluster, maxcluster, prevcluster, maxclusterid;
  int *clustersize;

  static int pixperline = 32; // PGM format limits to 70 characters per line

  if (ncluster > MAXNCLUSTER)
    {
      printf("percwrite: WARNING ncluster too large, resetting to %d\n",
	     MAXNCLUSTER);

      ncluster = MAXNCLUSTER;
    }

  if (ncluster > 1)
    {
      printf("percwrite: visualising the largest %d clusters\n", ncluster);
    }
  else
    {
      printf("percwrite: only visualising the largest cluster\n");
    }

  /*
   * Allocate the local clustersize array
   */

  maxclusterid = l*l;

  if ((clustersize = (int *) malloc((maxclusterid+1)*sizeof(int))) == NULL)
    {
      printf("percwrite: allocation of clustersize failed\n");
      exit(1);
    }

  /*
   * Count up the size of each cluster
   */

  for (i=0; i<=maxclusterid; i++)
    {
      clustersize[i] = 0;
    }

  for (i=0; i < l; i++)
    {
      for (j=0; j < l; j++)
	{
	  clusterid = map[i][j];

	  if (clusterid > 0)
	    {
	      clustersize[clusterid]++;
	    }
	}
    }

  /*
   * Find the size of the "ncluster" largest clusters (by brute force!)
   */

  prevcluster = maxclusterid+1; // Larger than the largest possible cluster id

  for (icluster = 0; icluster < ncluster; icluster++)
    {
      maxcluster = 0;

      for (i=0; i<=maxclusterid; i++)
	{
	  if (clustersize[i] > maxcluster && clustersize[i] < prevcluster)
	    {
	      maxcluster = clustersize[i];
	    }
	}

      foundcluster[icluster] = maxcluster;
      prevcluster = maxcluster;
    }

  if (ncluster > 1)
    {
      printf("percwrite: cluster sizes are ");
    }
  else
    {
      printf("percwrite: maximum cluster size is ");
    }

  for (icluster = 0; icluster < ncluster-1; icluster++)
    {
      printf("%d, ", foundcluster[icluster]);
    }
  printf("%d\n", foundcluster[ncluster-1]);

  /*
   *  Write the file
   */

  printf("percwrite: opening file <%s>\n", percfile);

  fp = fopen(percfile, "w");

  printf("percwrite: writing data ...\n");

  /*
   *  Start with the PGM header
   */

  fprintf(fp, "P2\n");
  fprintf(fp, "%d %d\n%d\n", l, l, ncluster);

  /*
   *  Now write the cells to file so that map[0][0] is in the
   *  bottom-left-hand corner and map[l-1][l-1] is in the
   *  top-right-hand corner
   */

  npix = 0;

  for (j=l-1; j >= 0; j--)
    {
      for (i=0; i < l; i++)
	{
	  clusterid = map[i][j];

	  /*
	   * Write out the largest cluster(s), shading appropriately
	   */

	  colour = 0;

	  if (clusterid > 0)
	    {
	      for (icluster = 0; icluster < ncluster; icluster++)
		{
		  if (clustersize[clusterid] == foundcluster[icluster])
		    {
		      // Largest (first) cluster is white

		      colour = ncluster - icluster;
		    }
		}
	    }

	  npix++;

	  // Make sure lines wrap after "npix" pixels

	  if (npix == 1)
	    {
	      fprintf(fp, "%1d", colour);
	    }
	  else if (npix < pixperline)
	    {
	      fprintf(fp, " %1d", colour);
	    }
	  else
	    {
	      fprintf(fp, " %1d\n", colour);
	      npix = 0;
	    }
	}
    }

  if (npix != 0) fprintf(fp, "\n");

  printf("percwrite: ... done\n");

  fclose(fp);
  printf("percwrite: file closed\n");

  /*
   * De-allocate the local clustersize array
   */

  free(clustersize);
}
