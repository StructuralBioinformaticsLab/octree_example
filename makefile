_processor = $(shell uname -p)

CC = icc
MOL_VERSION = 0.0.6
MOL_INCLUDE ="\"minilibmol/mol.$(MOL_VERSION).h\""
CPPFLAGS := -D _MOL_VERSION_="\"$(MOL_VERSION)\"" -D _MOL_INCLUDE_=$(MOL_INCLUDE) $(CPPFLAGS)
CPPFLAGS := -D ATOM_PRM="\"atom.$(MOL_VERSION).prm\"" $(CPPFLAGS)

LDFLAGS = -L$(HOME)/lib -L/opt/local/lib  -L/usr/local/lib/
CPPFLAGS := -I$(HOME)/include -I/opt/local/include $(CPPFLAGS)

CFLAGS =  -O3  -AVX -xhost -parallel -ip -vec-report -ansi-alias -restrict -DVECTORIZE_OCTREE -Wall -W -Wshadow -Wpointer-arith -Wcast-qual
LIBS := -Lminilibmol -lmol.$(MOL_VERSION) -lm -lpapi $(LIBS)

OBJS = test.o 

test:	$(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -Xlinker -zmuldefs -o test

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

all:
	test
test.mpi: $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(OBJS)
	$(RM) test
