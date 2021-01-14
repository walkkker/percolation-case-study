MF=	Makefile

CC=	mpicc
CFLAGS=	-cc=icc -O3 -Wall

LFLAGS= $(CFLAGS)

EXE=	percolate

INC= \
	percolate.h

SRC= \
	percolate.c \
	percio.c \
	unirand.c \
	test_percolation.c \
	swap_halos.c \
	distributing_map.c \
	collecting_map.c \
	cart_topology.c



#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(OBJ):	$(INC)

$(EXE):	$(OBJ)
	$(CC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ):	$(MF)

clean:
	rm -f $(EXE) $(OBJ) core
