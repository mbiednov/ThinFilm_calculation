#include "mex.h"
#include "multilayer.cpp"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *inMatrix_n;               /* real part of index 1xN input matrix */
    double *inMatrix_k;               /* imagnary part of 1xN input matrix */
    double *inMatrix_thickness;       /* thickness of 1xN input matrix */
    double *outMatrix;              /* output matrix */
    
    //mexPrintf("the input accepts as: \nint NumberLayers, array n, array k, array thickness, double angleIncidence(deg), double wavelength(nm), bool P(true)_S(false)\n");
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 7) {
        mexErrMsgIdAndTxt( "MATLAB:mxcreatecellmatrix:minrhs",
                "the input must be int NumberLayers, array n, array k, array thickness, double angleIncidence(deg), double wavelength(nm),bool P(true)_S(false)\n");
    } 
    if(nlhs > 1){
        mexErrMsgIdAndTxt( "MATLAB:mxcreatecellmatrix:maxlhs",
                "Too many output arguments.");
    }
    
    /* get the number of the layers  */
    int numofLayers = mxGetScalar(prhs[0]);
    
    /* create a pointer to the real data in the input matrix  */
    inMatrix_n = mxGetPr(prhs[1]);
    inMatrix_k = mxGetPr(prhs[2]);
    inMatrix_thickness = mxGetPr(prhs[3]);
    
    // size of n and k array
    int ncols_n = mxGetN(prhs[1]);
    int ncols_k = mxGetN(prhs[2]);
    int ncols_thickness = mxGetN(prhs[3]);
    
    //for (int i0=0;i0<ncols_n;i0++)
    //{
    //    mexPrintf("%.2f\t%.2f\t%.2f\t\n",inMatrix_n[i0],inMatrix_k[i0],inMatrix_thickness[i0]);
    //}
    
    //
    double angleofIncidence = mxGetScalar(prhs[4]);
    double wavelength = mxGetScalar(prhs[5]);
    bool SPindicator = mxGetScalar(prhs[6]);
    
    if (ncols_n!=ncols_k ||ncols_n!=ncols_thickness)
    {
        mexErrMsgIdAndTxt( "MATLAB:mxcreatecellmatrix:minrhs",
                "the input arraies are not equal.");
    }
    
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    outMatrix = mxGetPr(plhs[0]);
    
    outMatrix[0] = calculateReflectivity(numofLayers, inMatrix_n,inMatrix_k, inMatrix_thickness, angleofIncidence, wavelength, SPindicator);
    
    
    
}