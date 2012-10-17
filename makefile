_processor = $(shell uname -p)

MOL_VERSION = 0.0.6
MOL_INCLUDE ="\"mol.$(MOL_VERSION).h\""

CPPFLAGS := -D _MOL_VERSION_="\"$(MOL_VERSION)\"" -D _MOL_INCLUDE_=$(MOL_INCLUDE) $(CPPFLAGS)
CPPFLAGS := -D ATOM_PRM="\"atom.$(MOL_VERSION).prm\"" $(CPPFLAGS)

LDFLAGS := -L$(HOME)/lib -L/opt/local/lib  -L/usr/local/lib -Wl,-rpath,/usr/local/lib
CPPFLAGS := -I$(HOME)/include -I/opt/local/include -I/usr/local/include $(CPPFLAGS)

CFLAGS =  -O3 -Wall -W -Wshadow -Wpointer-arith -Wcast-qual
LIBS := -lmol.$(MOL_VERSION) -lm $(LIBS)

mintest: 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c mintest.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) mintest.o $(LIBS) -o mintest

papi: 
	$(CC) $(CPPFLAGS) -D USE_PAPI $(CFLAGS) -c mintest.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) mintest.o $(LIBS) -lpapi -o mintest

all:
	mintest

test: 
	./mintest --pdb examples/1acb_b_rmin.pdb --vdw
	./mintest --pdb examples/1acb_b_rmin.pdb --nofixed --vdw
	./mintest --pdb examples/1acb_b_rmin.pdb --nofixed --vdw --hbond	
	./mintest --pdb examples/1acb_b_rmin.pdb --nofixed --vdw --vdwcut 14 18 --hbond
	./mintest --pdb examples/1acb_b_rmin.pdb --nofixed --vdw --vdwaprx 0.5 --hbond

clean:
	$(RM) mintest.o
	$(RM) mintest
