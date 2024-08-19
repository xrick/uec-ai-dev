/*
 Mex-file by:
 ------- kbdwin.m -----------------------------------------
 Marios Athineos, marios@ee.columbia.edu
 http://www.ee.columbia.edu/~marios/
 Copyright (c) 2002 by Columbia University.
 All rights reserved.
 ----------------------------------------------------------

 Original C code by:
 Programmer:    Craig Stuart Sapp <craig@ccrma.stanford.edu>
 Creation Date: Sat Jan 27 14:27:14 PST 2001
 Last Modified: Sat Jan 27 15:24:58 PST 2001
 Filename:      kbdwindow.cpp
 $Smake:        g++ -O6 -o %b %f
 Syntax:        C++; functions in ANSI C

 Description:   This is a sample program for generating
                Kaiser-Bessel Derived Windows for audio signal
                processing -- in particular for the Time-Domain
                Alias Cancellation procedure which has the
                overlap-add requirements:
                   Window_m[N-1-n] + Window_m+1[n]^2 = 1.0;
                which means: The squares of the overlapped
                windows must add to the constant value 1.0
                for time domain alias cancellation to work.

                The two function necessary to create the KBD window are:
                   KBDWindow -- calculates the window values.
                   BesselI0  -- bessel function needed for KBDWindow.
*/

#include "mex.h"  /* for mex functions  */
#include <math.h> /* for pow function  */

/* function declarations */
void KBDWindow(double* window, int size, double alpha);
double BesselI0(double x);

/*
 Entry-point for Matlab calls
*/
void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    /* This is the syntax of our function */
    /* y = kbdwin(len,alpha) */
    double *y;
    double  len,alpha;
    
    /* Check arguments */
    if (nrhs != 2) { 
        mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments."); 
    } else if (!mxIsNumeric(prhs[0])) {
        mexErrMsgTxt("Argument must be numeric.");
    } else if (!mxIsNumeric(prhs[1])) {
        mexErrMsgTxt("Argument must be numeric.");
    } else if (mxGetNumberOfElements(prhs[0]) != 1 || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Argument must be non-complex scalar.");
    } else if (mxGetNumberOfElements(prhs[1]) != 1 || mxIsComplex(prhs[1])) {
        mexErrMsgTxt("Argument must be non-complex scalar.");
    }

    /* get the scalar value of the inputs len and alpha */
    len   = mxGetScalar(prhs[0]);
    alpha = mxGetScalar(prhs[1]);

    /* create a len-by-1 matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL);

    /* assign a pointer to the output */
    y = mxGetPr(plhs[0]);
    
    /* call the timestwo_alt subroutine */
    KBDWindow(y,len,alpha);
}

/*
 KBDWindow -- Kaiser Bessel Derived Window
      fills the input window array with size samples of the
      KBD window with the given tuning parameter alpha.
*/

#ifndef PI
#define PI 3.14159265358979323846264338328
#endif

void KBDWindow(double* window, int size, double alpha) {
    double sumvalue = 0.0;
    int i;
    
    for (i=0; i<size/2; i++) {
        sumvalue += BesselI0(PI * alpha * sqrt(1.0 - pow(4.0*i/size - 1.0, 2)));
        window[i] = sumvalue;
    }
    
    /* need to add one more value to the nomalization factor at size/2: */
    sumvalue += BesselI0(PI * alpha * sqrt(1.0 - pow(4.0*(size/2)/size-1.0, 2)));
    
    /* normalize the window and fill in the righthand side of the window: */
    for (i=0; i<size/2; i++) {
        window[i] = sqrt(window[i]/sumvalue);
        window[size-1-i] = window[i];
    }
}

/*
 BesselI0 -- Regular Modified Cylindrical Bessel Function (Bessel I).
*/
double BesselI0(double x) {
    double denominator;
    double numerator;
    double z;
    
    if (x == 0.0) {
        return 1.0;
    } else {
        z = x * x;
        numerator = (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* 
            (z* 0.210580722890567e-22  + 0.380715242345326e-19 ) +
            0.479440257548300e-16) + 0.435125971262668e-13 ) +
            0.300931127112960e-10) + 0.160224679395361e-7  ) +
            0.654858370096785e-5)  + 0.202591084143397e-2  ) +
            0.463076284721000e0)   + 0.754337328948189e2   ) +
            0.830792541809429e4)   + 0.571661130563785e6   ) +
            0.216415572361227e8)   + 0.356644482244025e9   ) +
            0.144048298227235e10);
        
        denominator = (z*(z*(z-0.307646912682801e4)+
            0.347626332405882e7)-0.144048298227235e10);
    }
    
    return -numerator/denominator;
}
