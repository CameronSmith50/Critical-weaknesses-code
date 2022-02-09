/**************************************************************************
 * index = quickfind(input,sum,randnum)
 *************************************************************************/

#include "mex.h"
#include "math.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *input, *index, *parameter, randnum, cumsum, sum;
    mwSize mrows, ncols, total, i;
    
    input = mxGetPr(prhs[0]); /* Pointer to the input matrix */
    parameter= mxGetPr(prhs[1]);
    sum= *parameter; /* sum of all values in input */
    parameter= mxGetPr(prhs[2]);
    randnum= *parameter; /* random number */
    randnum*=sum; /* scale by sum */
    
    /* Dimensions of input matrix */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    total = mrows*ncols;
    
    /* Create output */
    /*   plhs[0] = mxCreateScalarDouble(0); */ /* Create output */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    index = mxGetPr(plhs[0]);
    index[0] = 0;
    
    /* Find chosen entry of input */
    cumsum=0;
    for(i=0;i<total;i++){
        cumsum+=input[i];
        if(cumsum>=randnum){
            index[0]=i+1;
            break;
        }
    }
}
