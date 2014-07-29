octree_example
==============

Installation
------------
Requirements:
	- icc (or with modifications to the makefile, a version of gcc that has cilk support)
	- papi ( http://icl.cs.utk.edu/papi/ )

You should first build libmol, then this program

    cd minilibmol
    make
    cd ..
    make


Ths program require the CHARMM19 parameters as found at http://mackerell.umaryland.edu/CHARMM_ff_params.html

You can download the latest prms (e.g. toppar_c36_dec13.tgz ), extract the directory.

Then, run:

    cp toppar/toph19.inp params/pdbamino.rtf
    cat toppar/param19.inp toppar/ace/acepar19.inp > params/parm.prm
