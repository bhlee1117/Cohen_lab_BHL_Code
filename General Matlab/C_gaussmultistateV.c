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
 * unsigned char *states (N elements)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *means = mxGetPr(prhs[0]);
    double *vars = mxGetPr(prhs[1]);
    double *transitions = mxGetPr(prhs[2]);
    double *intens = mxGetPr(prhs[3]);
    const mwSize numstates = mxGetM(prhs[0]);
    const mwSize numbins = mxGetM(prhs[3]);

    double *logvarterms = mxCalloc(numstates,sizeof(double));
    double *logtransitions = mxCalloc(numstates*numstates,sizeof(double));
    double *delta = mxCalloc(numstates,sizeof(double));
    double *tempdelta = mxCalloc(numstates*numstates,sizeof(double));
    unsigned char *states;
    unsigned char *backtracker = mxCalloc(numstates*numbins,sizeof(unsigned char));
    double maxlikelihood;
    mwIndex i,j,t;

/*
    mexPrintf("%d bins\n",numbins);
*/
    
    plhs[1] = mxCreateNumericMatrix(numbins,1,mxUINT8_CLASS,0); /*states matrix*/
    
    states = mxGetData(plhs[1]);
    
    if(states == NULL)
        mexErrMsgTxt("Not enough memory to create states vector!");
    
    for(i = 0;i < numstates;i++)
    {
        logvarterms[i] = -log(2*M_PI*vars[i])/2;
        for(j = 0;j < numstates;j++)
            logtransitions[j+numstates*i] = log(transitions[j+numstates*i]);
        delta[i] = log(0);
    }
    delta[0] = (intens[0]-means[0])*(means[0]-intens[0])/(2*vars[0])+logvarterms[0];
        /* initialize with full probability in state 0 */

    for(t = 1; t < numbins; t++)
    {
        for(i = 0; i < numstates; i++)
            for(j = 0; j < numstates; j++)
                tempdelta[j+numstates*i] = delta[i]+logtransitions[j+numstates*i];
        for(j = 0;j < numstates;j++)
        {
            backtracker[j+numstates*t] = (unsigned char)j; /*favor staying in the same state in case of tie*/
            maxlikelihood = tempdelta[j+numstates*j];
            for(i = 0; i < j; i++)
            {
                if(tempdelta[j+numstates*i] > maxlikelihood)
                {
                    backtracker[j+numstates*t] = (unsigned char)i;
                    maxlikelihood = tempdelta[j+numstates*i];
                }
            }
            for(i = j+1; i < numstates; i++) /*split into two loops to avoid repeating the same state*/
            {
                if(tempdelta[j+numstates*i] > maxlikelihood)
                {
                    backtracker[j+numstates*t] = (unsigned char)i;
                    maxlikelihood = tempdelta[j+numstates*i];
                }
            }
            delta[j] = maxlikelihood+(intens[t]-means[j])*(means[j]-intens[t])/(2*vars[j])+logvarterms[j];
        }
    }
    states[numbins-1] = 0;
    maxlikelihood = delta[0];
    for(i = 1;i < numstates;i++)
    {
        if(delta[i] > maxlikelihood)
        {
            states[numbins-1] = (unsigned char)i;
            maxlikelihood = delta[i];
        }
    }
    plhs[0] = mxCreateDoubleScalar(maxlikelihood);
    for(t = numbins-1;t >= 1; t--)
        states[t-1] = (unsigned char)backtracker[(mwIndex)states[t]+numstates*t];
    mxFree(logvarterms);
    mxFree(logtransitions);
    mxFree(delta);
    mxFree(tempdelta);
    mxFree(backtracker);
}
