#include "mex.h"
#include<math.h>

/* ARGUMENTS:
 * double *means (Mx1 elements)
 * double *vars (Mx1 elements)
 * double *transitions (MxM elements)
 * double *intens (Nx1 elements)
 * 
 * RETURNS:
 * double loglikelihood
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *means = mxGetPr(prhs[0]);
    double *vars = mxGetPr(prhs[1]);
    double *transitions = mxGetPr(prhs[2]);
    double *intens = mxGetPr(prhs[3]);
    const mwSize numstates = mxGetM(prhs[0]);
    const mwSize numbins = mxGetM(prhs[3]);

    double *probs = mxCalloc(numstates,sizeof(double));
    double *tempprobs = mxCalloc(numstates,sizeof(double));

    double loglikelihood = (intens[0]-means[0])*(means[0]-intens[0])/(2*vars[0])-log(2*M_PI*vars[0])/2;
        /* can put in explicit value here since we always start in state 0 */
    double templikelihood;

    mwIndex i,j,t;

    probs[0] = 1; /* initialize with full probability in state 0 */
    for(i = 1; i < numstates; i++)
        probs[i] = 0;
    for(t = 1; t < numbins; t++)
    {
        templikelihood = 0;
        for(j = 0; j < numstates; j++)
        {
            tempprobs[j] = 0;
            for(i = 0; i < numstates; i++)
            {
                tempprobs[j] += probs[i]*transitions[j+numstates*i];
            }
            tempprobs[j] *= exp((intens[t]-means[j])*(means[j]-intens[t])/(2*vars[j]))/sqrt(2*M_PI*vars[j]);
            templikelihood += tempprobs[j];
        }
        for(j = 0; j < numstates; j++)
            probs[j] = tempprobs[j]/templikelihood;
        loglikelihood += log(templikelihood);
    }
    mxFree(tempprobs);
    mxFree(probs);
    plhs[0] = mxCreateDoubleScalar(loglikelihood);
}
