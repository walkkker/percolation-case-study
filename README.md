# 2D decomposition parallel program for percolation

## Code structure

1. percolate.h : Header file.

2. percolate.c : The main C file. This file includes the main function of this program.

3. cart_topology.c : This file has only one function.
```
  * cart_topology(): Create a Cartesian topology for all processes and get the number of processes in each dimension.
```

4. distributing_map.c : This file has two functions.
```
  * distributing_map_to_old(): Using non-blocking communications to transfer data directly from map to the internal region of old.
  * zero_halos(): Zero the top, bottom, left and right halos of old array.
```

5. swap_halos.c : This file has only one function.
```
  * swap_halos(): Non-blocking communications with four neighbors to swap halos.
```

6. collecting_map.c : The file has only one function.
```
  * collecting_map_from_old(): Using non-blocking communications to collecte data from old array on each processor to the map array on controller.
```

7. test_percolation.c: The file has only one function.
```
  * test_percolation(): Test to see if percolation occurred.
```

8. unirand.c : Random number generator.
9. percio.c : Function to write a percolation map in greyscale Portable Grey Map (PGM) format.


## Build the code on Cirrus
The code must build with Intel compilers and the mpt/2.22 and intel-compilers-19 modules.
On Cirrus, run:

```console
$ module load intel-compilers-19

$ module load mpt/2.22
```

## Compilation

The program is builded completely via make. Makefile is included in this folder. To build the program run the following:

```console
$ make
```

If you want to delete all generated files, run:

```console
$ make clean
```

## Run

To run the Percolate program:


Run the percolate executable file on the login node of Cirrus:

```console
$ mpirun -n [NUMBER OF PROCESSES] ./percolate [SEED]
```
(0 <= seed <= 900,000,000)

For example, to run with 16 processes and a seed 6543:

```console
$ mpirun -n 16 ./percolate 6543
```

By default, the code runs with L = 432 x 432 and rho = 0.411.


### Changing problem size

If you want to change the problem size, you need to modify the value of L in the header file "percolate.h", where on line 9

```c
#define L 432
```


If you want to change the value of rho, you need to modify the source file "percolate.c", where on line 128

```c
rho = 0.411;
```


### PGM output files

A Portable Grey Map (PGM) file is output.

This file does not include the halo as the use of a halo is an implementation detail.

This file is plain-text so you can view it as you would any plain-text file e.g.:

```console
$ cat map.pgm
```

PGM files can be viewed graphically using ImageMagick commands as follows.

Cirrus users will need first need to run:

```console
$ module load ImageMagick
```

To view a PGM file, run:

```console
$ display map.pgm
```



## Restrictions on the running of the code

The number of processes, P, must evenly divide the grid size, L.

If P cannot evenly divide L, the program will output as follows and terminate:

```
percolate: ERROR L is not an exact multiple of the number of processes
```


