octree_example
==============

Installation
------------

You should first build libmol, then this program

    cd libmol
    make && make install
    cd ..
    make


Ths program require the CHARMM19 parameters as found at http://mackerell.umaryland.edu/CHARMM_ff_params.html

You can download the latest prms (e.g. toppar_c36_aug12.tgz ), extract the directory.

Then, run:
    cp toppar/toph19.inp params/pdbamino.rtf
    cat toppar/param19.inp toppar/ace/acepar19.inp > params/parm.prm
