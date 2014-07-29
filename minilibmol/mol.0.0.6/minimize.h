/*
Copyright (c) 2009-2012, Structural Bioinformatics Laboratory, Boston University
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
- Neither the name of the author nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific
  prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef _MOL_MINIMIZE_H_
#define _MOL_MINIMIZE_H_

/** \file minimize.h
        This file contains  functions
        for local structure minimization
*/

typedef struct {
    int     n;
    int     m;
    int     niter;       /* number of iterations so far                        */
    int     nfuns;       /* number of function evaluations so far              */
    int     iflag;
    int     diagco;
    int     iprint[2];   /* see the comment in lbfgs.f for usage of this field */
    double  eps;
    double  xtol;
    double *diag;
    double *w;
} lbfgs_t;

typedef enum
{
    MOL_LBFGS,
    MOL_CONJUGATE_GRADIENTS,
    MOL_POWELL,
} mol_min_method;


double max (double x, double y);

lbfgs_t* lbfgs_create(int n, int m, double eps);

int lbfgs_run(lbfgs_t* obj, double* x, double f, double* g);

void lbfgs_destroy(lbfgs_t* obj);

void linestep ( double* stx , double* fx , double* dx , double* sty , double* fy , double* dy , double* stp , double fp , double dp , short* brackt , double stpmin , double stpmax , int* info );

double ddot ( int n, double* dx, int ix0, int incx, double* dy, int iy0, int incy );

void lb1 ( int* iprint , int* iter , int* nfun , double gnorm , int n , int m , double* x , double f , double* g , double* stp , short finish );

void daxpy ( int n , double da , double* dx , int ix0, int incx , double* dy , int iy0, int incy);

//LBFGS minimizer for arbitrary function 
void lbfgsMin(double* orig, unsigned int maxIt, double tol, int ndim,  void* prms,void (*egfun)(int , double* , void* , double* , double*),double* min, double* fmim );
//Wrapper for atom group minimizers;
void minimize_ag(mol_min_method min_type, unsigned int maxIt, double tol, struct atomgrp* ag,void* minprms, void (*egfun)(int , double* , void* , double* , double*));

// Powell/dirPowell below

void bracket(double* orig, double* dir, double step,
             int ndim, void* prms,
             void (*egfun)(int , double* , void* , double* , double*),
             double *fb,  double *la, double *lb, double *lc);

void brent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double tol, unsigned int maxtimes,
           double*  min, double* fmim);

void dirbrent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double tol, int maxtimes,
           double*  min, double* fmim, double*  grad);

void powell(double* orig, double* directions, unsigned int maxIt, double tol,
            int ndim, void* prms,
            void (*egfun)(int , double* , void* , double*, double* ),
            double* min, double* fmim);

void dirMin(double* orig, unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim);

void limin(double* orig, double* dir, unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim);

void old_brent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double*, double* ),
           double tol, int maxtimes,
           double*  min, double* fmim);
// Powell/dirPowell above

#endif
