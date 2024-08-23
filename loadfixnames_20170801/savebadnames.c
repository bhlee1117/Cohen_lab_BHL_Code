/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Function:    savebadnames
 * Filename:    savebadnames.c
 * Programmer:  James Tursa
 * Version:     1.00
 * Date:        June 18, 2013
 * Copyright:   (c) 2013 by James Tursa, All Rights Reserved
 *
 *  This code uses the BSD License:
 *
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are 
 *  met:
 *
 *     * Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer in 
 *       the documentation and/or other materials provided with the distribution
 *      
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Building:
 *
 * SAVEBADNAMES is typically self building. That is, the first time you call it,
 * the savebadnames.m file recognizes that the mex routine needs to be compiled and
 * then the compilation will happen automatically.
 *
 * The usage is as follows (arguments in brackets [ ] are optional):
 *
 * Syntax
 *
 * savebadnames(FILENAME [,verbose])
 *
 *     FILENAME = The mat file to be saved
 *     verbose  = 1 or 0 (optional, causes name change log to be displayed)
 *
 * Description
 *
 * SAVEBADNAMES saves a mat file using intentionally invalid names.
 * This is a program to create a test mat file for the LOADFIXNAMES function.
 *
 * Change Log:
 * 2013/Jun/18 --> 1.00, Initial Release
 *
 ****************************************************************************/

// Includes -----------------------------------------------------------

#include "mex.h"
#include "mat.h"
#include <ctype.h>
#include <string.h>

// Macros -------------------------------------------------------------

#ifndef MWSIZE_MAX
#define  mwIndex        int
#define  mwSignedIndex  int
#define  mwSize         int
#endif

// Prototypes ---------------------------------------------------------

char *appendmat(char *, mwSize);

// Gateway Function ---------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray *mx;
	char *filename, *matfilename;
	char mode[] = "w";
	char *names[] = {"Name with spaces",
		             "Name with dollar sign $",
					 "Lots of bad stuff $^%*&",
		             "123 started with number",
					 "^(*& started with bad chars",
					 "123^*&( another bad one"};
	char matlabmat[] = "matlab.mat";
	MATFile *mfp;
	int verbose = 1;
	mwSize i, n;

/* Check arguments */

	if( nlhs > 0 ) {
		mexErrMsgTxt("Too many outputs.");
	}
	if( nrhs > 2 ) {
		mexErrMsgTxt("Too many inputs.");
	}
	if( !mxIsChar(prhs[0]) ) {
		mexErrMsgTxt("1st argument must be filename char string.");
	}
	if( mxGetNumberOfDimensions(prhs[0]) != 2 || (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1) ) {
		mexErrMsgTxt("1st argument must be filename char string.");
	}
	if( nrhs >= 2 ) {
		if( mxIsNumeric(prhs[1]) ) {
			verbose = (mxGetScalar(prhs[1]) != 0.0);
		} else {
			mexErrMsgTxt("Invalid 2nd argument.");
		}
	}
	if( nrhs == 0 ) {
		matfilename = matlabmat; /* If no filename given, use "matlab.mat" */
	} else {
		filename = mxArrayToString(prhs[0]); /* Get the filename as a C-string */
		n = mxGetNumberOfElements(prhs[0]);
		matfilename = appendmat(filename, n); /* Append ".mat" if necessary */
	}
	if( verbose ) {
		mexPrintf("Attempting to open mat file: %s\n",matfilename);
	}
	mfp = matOpen(matfilename, mode); /* Open the mat file */
	if( mfp == NULL ) {
		mexPrintf("Tried to open file %s\n",matfilename);
		if( matfilename != matlabmat ) mxFree(matfilename);
		mexErrMsgTxt("Unable to open file.");
	}
	if( matfilename != matlabmat ) mxFree(matfilename); /* Done with temp filename C-string */
	n = sizeof(names) / sizeof(char *);
	if( verbose ) {
	    mexPrintf("About to put %d variables ...\n",n);
	}
	mx = mxCreateDoubleMatrix( 2, 2, mxREAL );
	for( i=0; i<n; i++ ) {
		if( verbose ) {
			mexPrintf("Putting variable with bad name: <%s>\n",names[i]);
		}
		matPutVariable( mfp, names[i], mx ); /* Put a bunch of bad named variables into the file */
	}
	mxDestroyArray(mx);
	matClose(mfp);
}

/**************************************************************************************************/

char *appendmat(char *filename, mwSize n)
{
	char *fname = filename;
	char *mname, *matfilename;

	if( n < 5 || (strcmp(filename+n-4,".mat") != 0 && strcmp(filename+n-4,".MAT") != 0 && strcmp(filename+n-4,".Mat") != 0) ) {
		mname = matfilename = mxMalloc(n+5);
		while( *filename ) {
			*matfilename++ = *filename++;
		}
		*matfilename++ = '.';
		*matfilename++ = 'm';
		*matfilename++ = 'a';
		*matfilename++ = 't';
		*matfilename = '\0';
		mxFree(fname);
	} else {
		mname = filename;
	}
	return mname;
}
