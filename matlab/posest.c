/* ////////////////////////////////////////////////////////////////////////////////
// 
//  Matlab MEX file for posest
//  Copyright (C) 2015  Manolis Lourakis (lourakis **at** ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//////////////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <posest.h>

#include "mex.h"

/**
#define DEBUG
**/

#define MAX(A, B)     ((A)>=(B)? (A) : (B))


/* display printf-style error messages in matlab */
static void matlabFmtdErrMsgTxt(char *fmt, ...)
{
char  buf[256];
va_list args;

	va_start(args, fmt);
	vsprintf(buf, fmt, args);
	va_end(args);

  mexErrMsgTxt(buf);
}

/* display printf-style warning messages in matlab */
static void matlabFmtdWarnMsgTxt(char *fmt, ...)
{
char  buf[256];
va_list args;

	va_start(args, fmt);
	vsprintf(buf, fmt, args);
	va_end(args);

  mexWarnMsgTxt(buf);
}

/* matlab matrices are in column-major, this routine converts them to row major for posest */
static double *getTranspose(mxArray *Am)
{
int m, n;
double *At, *A;
register int i, j;

  m=mxGetM(Am);
  n=mxGetN(Am);
  A=mxGetPr(Am);
  At=mxMalloc(m*n*sizeof(double));

  for(i=0; i<m; i++)
    for(j=0; j<n; j++)
      At[i*n+j]=A[i+j*m];
  
  return At;
}

/*
 [rt, idxOutliers, foc]=posest(pts2D, pts3D, inlPcent, K, NLrefine, verbose);
        pts2D, pts3D are the matched 2D-3D point coordinates
        inlPcent is the expected percentage of inliers (>=0.5)
        K is the 3x3 camera intrinsic calibration matrix or a 1x2 vector with the image width and height
        inlPcent is the expected fraction of inliers in the input pairs
        NLrefine controls non-linear refinement, can be one of: 'norefine', 'repr_err', 'repr_err_mlsl', 'objspc_err'
        or empty (implying 'norefine', default)
        verbose is optional
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *Prhs[])
{
register int i;
int NLrefine, verbose, nparms=NUM_RTPARAMS, ret, nbOutliers, *idxOutliers;
double (*pts2D)[2], (*pts3D)[3], inlPcent, K[9], rtf[NUM_RTFPARAMS];
register double *pdbl;
mxArray **prhs=(mxArray **)&Prhs[0];
int nmatches, len, status;

  /* parse input args; start by checking their number */
  if(nrhs<4 && nrhs>6)
    matlabFmtdErrMsgTxt("posest: between 4 and 6 input arguments required (got %d).", nrhs);
  if(nlhs>3)
    matlabFmtdErrMsgTxt("posest: at most 3 output arguments returned (got %d).", nlhs);
    
  /** pts2D **/
  /* first argument must be a two-column matrix */
  if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetN(prhs[0])!=2)
    matlabFmtdErrMsgTxt("posest: first argument must be a two-column matrix (got %dx%d).", mxGetM(prhs[0]), mxGetN(prhs[0]));
  pts2D=(double (*)[2])getTranspose(prhs[0]);
  nmatches=mxGetM(prhs[0]);

  /** pts3D **/
  /* second argument must be a three-column matrix */
  if(!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1])!=3)
    matlabFmtdErrMsgTxt("posest: second argument must be a three-column matrix (got %dx%d).", mxGetM(prhs[1]), mxGetN(prhs[1]));
  if(mxGetM(prhs[1])!=nmatches)
    matlabFmtdErrMsgTxt("posest: two first arguments should have the same number of rows (got %d and %d).", nmatches, mxGetM(prhs[1]));
  pts3D=(double (*)[3])getTranspose(prhs[1]);

  /** inlPcent **/
  /* the third argument must be a scalar */
  if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
    mexErrMsgTxt("posest: inlPcent must be a scalar.");
  inlPcent=mxGetScalar(prhs[2]);

  /** K **/
  /* check whether the fourth argument is a 3x3 matrix or 1x2 vector */
  if(mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3]) && mxGetM(prhs[3])==3 && mxGetN(prhs[3])==3){ // matrix
    pdbl=mxGetPr(prhs[3]);
    /* get transposed, i.e. row major */
    K[0]=pdbl[0]; K[1]=pdbl[3]; K[2]=pdbl[6];
    K[3]=pdbl[1]; K[4]=pdbl[4]; K[5]=pdbl[7];
    K[6]=pdbl[2]; K[7]=pdbl[5]; K[8]=pdbl[8];
  }
  else if(mxIsDouble(prhs[3]) && !mxIsComplex(prhs[3]) && mxGetM(prhs[3])==1 && mxGetN(prhs[3])==2){ // vector
    nparms=NUM_RTFPARAMS;
    pdbl=mxGetPr(prhs[3]);
    /* create a K matrix with zeroes for the focal lengths and principal point on the image center */
    K[0]=0.0; K[1]=0.0; K[2]=pdbl[0]*0.5;
    K[3]=0.0; K[4]=0.0; K[5]=pdbl[1]*0.5;
    K[6]=0.0; K[7]=0.0; K[8]=1.0;
  }
  else
    matlabFmtdErrMsgTxt("posest: fourth arguments should be either a 3x3 matrix or 1x2 vector (got %d x %d).", mxGetM(prhs[3]), mxGetN(prhs[3]));

  /** NLrefine **/
  /* check whether fifth argument is a string */
  if(nrhs>=5 && mxIsChar(prhs[4])==1 && mxGetM(prhs[4])==1){
    char *str;

    /* examine supplied name */
    len=mxGetN(prhs[4])+1;
    str=mxCalloc(len, sizeof(char));
    status=mxGetString(prhs[4], str, len);
    if(status!=0)
      mexErrMsgTxt("posest: not enough space. String is truncated.");

    for(i=0; str[i]; ++i)
      str[i]=tolower(str[i]);

    if(!strcmp(str, "norefine")) NLrefine=POSEST_REPR_ERR_NO_NLN_REFINE;
    else if(!strcmp(str, "repr_err")) NLrefine=POSEST_REPR_ERR_NLN_REFINE;
    else if(!strcmp(str, "repr_err_mlsl")) NLrefine=POSEST_REPR_ERR_NLN_MLSL_REFINE;
    else if(!strcmp(str, "objspc_err")) NLrefine=POSEST_OBJSPC_ERR_LHM;
    else matlabFmtdErrMsgTxt("posest: unknown minimization type '%s'.", str);

    if(NLrefine==POSEST_OBJSPC_ERR_LHM && nparms==NUM_RTFPARAMS)
      matlabFmtdErrMsgTxt("posest: focal length cannot be estimated in combination with 'objspc_err'.");

    mxFree(str);

    ++prhs;
    --nrhs;
  }
  else
    NLrefine=POSEST_REPR_ERR_NO_NLN_REFINE;

  /** verbose **/
  /* the sixth argument must be a scalar */
  if(nrhs>=5){
    if(!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetM(prhs[4])!=1 || mxGetN(prhs[4])!=1)
      mexErrMsgTxt("posest: verbose must be a scalar.");
    verbose=mxGetScalar(prhs[4])!=0.0;
  }
  else
    verbose=0;

  if(nlhs>1) /* outlier indices should be returned */
    idxOutliers=mxMalloc(nmatches*sizeof(int));
  else
    idxOutliers=NULL;

  /* invoke posest */
  ret=posest(pts2D, pts3D, nmatches, inlPcent, K, rtf, nparms, NLrefine, idxOutliers, &nbOutliers, verbose);
  if(ret!=POSEST_OK)
    mexErrMsgTxt("posest: did not complete successfully.");

  /* copy back returned results */
  /** rt **/
  plhs[0]=mxCreateDoubleMatrix(NUM_RTPARAMS, 1, mxREAL);
  pdbl=mxGetPr(plhs[0]);
  pdbl[0]=rtf[0]; pdbl[1]=rtf[1]; pdbl[2]=rtf[2];
  pdbl[3]=rtf[3]; pdbl[4]=rtf[4]; pdbl[5]=rtf[5];

  /** idxOutliers **/
  if(nlhs>1){
    plhs[1]=mxCreateDoubleMatrix(nbOutliers, 1, mxREAL);
    pdbl=mxGetPr(plhs[1]);
    for(i=0; i<nbOutliers; ++i)
      *pdbl++=idxOutliers[i]+1; // convert from 0 to 1-based

    mxFree(idxOutliers);
  }

  /** focal **/
  if(nlhs>2){
    plhs[2]=mxCreateDoubleMatrix(1, 1, mxREAL);
    pdbl=mxGetPr(plhs[2]);
    if(nparms==NUM_RTFPARAMS){ // estimated f
      pdbl[0]=rtf[6];
    }
    else
      pdbl[0]=K[0]; // f from supplied intrinsics
  }

  /* cleanup */
  mxFree(pts2D);
  mxFree(pts3D);
}
