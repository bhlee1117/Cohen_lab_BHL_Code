/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Function:    loadfixnames
 * Filename:    loadfixnames.c
 * Programmer:  James Tursa
 * Version:     1.20
 * Date:        October 31, 2013
 * Copyright:   (c) 2013, 2017 by James Tursa, All Rights Reserved
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
 * LOADFIXNAMES is typically self building. That is, the first time you call it,
 * the loadfixnames.m file recognizes that the mex routine needs to be compiled and
 * then the compilation will happen automatically.
 *
 * The usage is as follows (arguments in brackets [ ] are optional):
 *
 * Syntax
 *
 * loadfixnames(FILENAME [,names] [,verbose])
 * S = loadfixnames(FILENAME [,names] [,verbose])
 *
 *     FILENAME = The mat file to be loaded
 *     names    = string or cell array of strings (variable name(s) to load)
 *     verbose  = 1 (optional, causes name change log to be displayed)
 *     S = struct returned instead of loading variables into workspace
 *
 * Description
 *
 * LOADFIXNAMES loads a mat file into the workspace, fixing invalid names.
 * All invalid characters are replaced with an underscore. Also, if the
 * first character is not a letter, an 'A' is added to the front. If the
 * variable is a structure, fixes the field names also. In the case of
 * field names, the name length is kept constant. So if a field name
 * begins with a digit it will be replaced with 'A' - 'J' instead. If the
 * field name begins with an invalid non-digit it will be replaced with
 * 'A' - 'Z' or 'a' - 'z' (letters are cycled in an attempt to avoid
 * name clashes). If the mat file contains duplicate variable names,
 * then LOADFIXNAMES will append a 5-digit number to the end of the
 * duplicate names, starting with 00002 and increasing by 1 each time.
 *
 * Limitations: The current renaming scheme makes a mild attempt to avoid
 *              name clashes, but does not guarantee this.
 *
 * See the companion function SAVEBADNAMES for a test function that will
 * intentionally create a mat file with invalid names.
 *
 * Change Log:
 * 2013/Jun/06 --> 1.00, Initial Release to Answers (used deep copy)
 * 2013/Jun/18 --> 1.10, Initial Release to FEX (uses shared data copy)
 * 2013/Oct/31 --> 1.20, Added code to deal with duplicate variable names
 *                       Added capability to return a struct
 * 2017/Aug/01 --> 1.30, Fixed a bug in a name string copy
 *                       Added smarter replacement characters
 *
 ****************************************************************************/

/* Includes ----------------------------------------------------------- */

#include "mex.h"
#include "mat.h"
#include <ctype.h>
#include <string.h>

/* Macros ------------------------------------------------------------- */

#if !defined(MWSIZE_MAX)
#if !defined(mwIndex)
#define mwIndex int
#endif
#if !defined(mwSize)
#define mwSize int
#endif
#endif

#define  MAXNAME  63

/* Prototypes --------------------------------------------------------- */

void fixvarname(char *, char *, int);
void fixfieldnames(mxArray *, int);
char *appendmat(char *, mwSize);
int isamatch(const char *, mxArray *);
int wildcmp(const char *wild, const char *str);

/* Global variables --------------------------------------------------- */

char alphabet[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
char alphanum[] = "AaABbBCcCDdDEeEFfFGgGHhHIiIJjJKkKLlLMmMNnNOoOPpPQqQRrRSsSTtTUuUVvVWwWXxXYyYZzZ_0a1b2c3d4e5f6g7h8i9j";
char upperASCII[] = "E__f__TT_PS_O_Z__________TS_o_zY_icLcYIS_C____R_o_23_uP__1o_423QAAAAAAACEEEEIIIIDNOOOOOx0UUUUYbBaaaaaaaceeeeiiiionooooo_ouuuuyby";

/* Prototype for undocumented API function --------------------------------- */

#ifndef MXCREATESHAREDDATACOPY_H
#define MXCREATESHAREDDATACOPY_H
mxArray *mxCreateSharedDataCopy(const mxArray *);
#endif

/* mxArray Header Structure definition ------------------------------------- */

#ifndef MXARRAYSTRUCT_H
#define MXARRAYSTRUCT_H

struct mxArrayStruct {
    void *name;             /*   prev - R2008b: Name of variable in workspace
				               R2009a - R2010b: NULL
				               R2011a - later : Reverse CrossLink pointer    */
    mxClassID ClassID;      /*  0 = unknown     10 = int16
                                1 = cell        11 = uint16
                                2 = struct      12 = int32
                                3 = logical     13 = uint32
                                4 = char        14 = int64
                                5 = void        15 = uint64
                                6 = double      16 = function_handle
                                7 = single      17 = opaque (classdef)
                                8 = int8        18 = object (old style)
                                9 = uint8       19 = index (deprecated)
                               10 = int16       20 = sparse (deprecated)     */
    int VariableType;       /*  0 = normal
                                1 = persistent
                                2 = global
                                3 = sub-element (field or cell)
                                4 = temporary
                                5 = (unknown)
                                6 = property of opaque class object
                                7 = (unknown)                                */
    mxArray *CrossLink;     /* Address of next shared-data variable          */
    size_t ndim;            /* Number of dimensions                          */
    unsigned int RefCount;  /* Number of extra sub-element copies            */
    unsigned int flags;     /*  bit  0 = is scalar double full
                                bit  2 = is empty double full
                                bit  4 = is temporary
                                bit  5 = is sparse
                                bit  9 = is numeric
                                bits 24 - 31 = User Bits                     */
    union {
        size_t M;           /* Row size for 2D matrices, or                  */
        size_t *dims;       /* Pointer to dims array for nD > 2 arrays       */
    } Mdims;
    size_t N;               /* Product of dims 2:end                         */
    void *pr;               /* Real Data Pointer (or cell/field elements)    */
    void *pi;               /* Imag Data Pointer (or field information)      */
    union {
        mwIndex *ir;        /* Pointer to row values for sparse arrays       */
        char *ClassName;    /* Pointer to Old User Defined Class Name        */
        mxClassID ClassID;  /* New User Defined Class ID (classdef)          */
    } irClassNameID;
    union {
        mwIndex *jc;        /* Pointer to column values for sparse arrays    */
        mxClassID ClassID;  /* Old User Defined Class ID                     */
    } jcClassID;
    size_t nzmax;           /* Number of elements allocated for sparse       */
/*  size_t reserved;           Don't believe this! It is not really there!   */
};

#endif

/* -------------------------------------------------------------------- */

int mexPutVariableSharedDataCopy(const char *workspace, const char *varname, const mxArray *value)
{
	void *vp;
	char *cp;
	mxArray *small, *value_sdc, *my;
	const mxArray *PropertyPtr;
	struct mxArrayStruct *mpv, *mpp;
	struct mxArrayStruct tpv;
	int double_linked_list, name_string, ipv, ipp, f;


	/* Check for non-NULL inputs */
	if( !workspace || !varname || !value ) return 1;

	/* Check for proper variable name */
	if( !isalpha(*varname) ) return 1; /* First character must be alpha */
	cp = varname;
	f = MAXNAME; /* Max variable name length in MATLAB */
	while( *++cp ) {
		if( !--f || !(isalnum(*cp) || *cp == '_') ) return 1; /* Remaining chars must be alphanumeric or underscore */
	}

	/* Put a small variable into the workspace with desired name */
	small = mxCreateDoubleMatrix( 1, 1, mxREAL );
	f = mexPutVariable( workspace, varname, small );
	mxDestroyArray( small );
	if( f ) return 1;

/* Remaining code is mostly borrowed from mxSetPropertySharedDataCopy */

	/* Get the pointer to the workspace variable */
	PropertyPtr = mexGetVariablePtr( workspace, varname );
	if( PropertyPtr == NULL ) return 1;
	mpp = (struct mxArrayStruct *) PropertyPtr;

	/* Make a shared data copy of the input mxArray */
	value_sdc = mxCreateSharedDataCopy( value );
	mpv = (struct mxArrayStruct *) value_sdc;

	/* Determine if double linked list is present (R2011a and later) */
	double_linked_list = mpv->name == value;

	/* Determine if first pointer is a name pointer (thru R2008b) */
	name_string = !double_linked_list && mpp->name != NULL;

	/* Swapping contents of PropertyPtr with value_sdc */
	tpv = *mpv; *mpv = *mpp; *mpp = tpv;

	/* Swap back the variable type, which we do not want swapped */
	ipv = mpv->VariableType;
	ipp = mpp->VariableType;
	mpv->VariableType = ipp;
	mpp->VariableType = ipv;

	/* Swap back the variable name pointer if appropriate */
	if( name_string ) {
		vp = mpp->name; mpp->name = mpv->name; mpv->name = vp;
	}

	/* Fix the flags, resetting the temporary flag */
	mpp->flags &= ~0x10u;

	/* Fix the CrossLink list so workspace variable is in the list correctly */
	my = mpp->CrossLink;
	while( my ) {
		mpp = (struct mxArrayStruct *) my;
		if( double_linked_list ) { /* Fix the reverse CrossLink pointer */
			double_linked_list = 0;
			mpp->name = PropertyPtr;
		}
		if( (my = mpp->CrossLink) == value_sdc ) {
			mpp->CrossLink = PropertyPtr;
			break;
		}
	}
	mxDestroyArray( value_sdc ); /* Our small variable, not the shared data copy anymore */
	return 0;
}

// Gateway Function ---------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray **mxarrays;
	mxArray *mx, *compare = NULL;
	char *filename, *matfilename;
	char mode[] = "r";
	const char *name;
	char *cname, *source, *target;
	char newname[MAXNAME+2];
	char matlabmat[] = "matlab.mat";
	MATFile *mfp;
	int e, j, k, duplicate, verbose = 0;
	mwSize n, ndim;
	mwSize *dims;
    char **names;
    int knames, nnames;

/* Check arguments */

	if( nlhs > 1 ) {
		mexErrMsgTxt("Too many outputs.");
	}
	if( nrhs > 3 ) {
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
		} else if( mxIsChar(prhs[1]) || mxIsCell(prhs[1]) ) {
			compare = prhs[1];
		} else {
			mexErrMsgTxt("Invalid 2nd argument.");
		}
	}
	if( nrhs >= 3 ) {
		if( mxIsNumeric(prhs[2]) ) {
			verbose = (mxGetScalar(prhs[2]) != 0.0);
		} else if( mxIsChar(prhs[2]) || mxIsCell(prhs[2]) ) {
			compare = prhs[2];
		} else {
			mexErrMsgTxt("Invalid 3rd argument.");
		}
	}
	if( verbose ) {
		mexPrintf("\nVerbose ...\n\n");
	}
	if( nrhs == 0 ) {
		matfilename = matlabmat; /* If no filename given, use "matlab.mat" */
	} else {
		filename = mxArrayToString(prhs[0]); /* Get the filename as a C-string */
		n = mxGetNumberOfElements(prhs[0]);
		matfilename = appendmat(filename, n); /* Append ".mat" if necessary */
	}
	mfp = matOpen(matfilename, mode); /* Open the mat file */
	if( mfp == NULL ) {
		mexPrintf("Tried to open file %s\n",matfilename);
		if( matfilename != matlabmat ) mxFree(matfilename);
		mexErrMsgTxt("Unable to open file.");
	}
	if( matfilename != matlabmat ) mxFree(matfilename); /* Done with temp filename C-string */
	mxarrays = NULL;
	names = NULL;
	knames = 0;
	nnames = 0;
	while( 1 ) {
		mx = matGetNextVariable( mfp, &name );
		if( mx == NULL ) break;
		if( !isamatch(name,compare) ) continue;
		fixvarname( name, newname, verbose );
		fixfieldnames( mx, verbose );
		if( knames == nnames ) {
			names = mxRealloc(names,(nnames+=1024)*sizeof(*names));
			if( nlhs ) {
			    mxarrays = mxRealloc(mxarrays,(nnames+=1024)*sizeof(*mxarrays));
			}
		}
		names[knames] = mxMalloc(MAXNAME+2+5);
		source = newname;
		target = names[knames];
		while( *target++ = *source++ ) {  /* Fixed 1.30 */
		}
		e = source - newname;
		duplicate = 0;
DUPLICATE:
		for( k=0; k<knames; k++ ) {
			if( strcmp(names[k],names[knames]) == 0 ) {
				if( duplicate ) {
					j = e + 4;
					while( names[knames][j] == '9' ) {
						names[knames][j--] = '0';
						if( j < e ) {
			                mexWarnMsgTxt("Too many duplicate names");
							goto TOOMANY;
						}
					}
					names[knames][j]++;
				} else {
					names[knames][e  ] = '0';
					names[knames][e+1] = '0';
					names[knames][e+2] = '0';
					names[knames][e+3] = '0';
					names[knames][e+4] = '2';
					names[knames][e+5] = '\0';
				    duplicate = 1;
				}
				goto DUPLICATE;
			}
		}
TOOMANY:
		if( duplicate && verbose ) {
			mexPrintf("Duplicate name <%s> renamed to <%s>\n",newname,names[knames]);
		}
		if( verbose ) {
			mexPrintf("Loading variable <%s> %s ",names[knames],mxGetClassName(mx));
			ndim = mxGetNumberOfDimensions(mx);
			dims = mxGetDimensions(mx);
			for( j=0; j<ndim-1; j++ ) {
				mexPrintf("%d x ",dims[j]);
			}
			mexPrintf("%d\n",dims[ndim-1]);
		}
		if( nlhs ) {
			mxarrays[knames] = mx;
		} else {
		    if( mexPutVariableSharedDataCopy( "caller", names[knames], mx ) ) {
			    mexPrintf("Tried to put variable into caller workspace: %s\n",names[knames]);
			    mexWarnMsgTxt("Unable to put variable into caller workspace.");
		    }
		    mxDestroyArray(mx);
		}
		knames++;
	}
	matClose(mfp);
	if( nlhs ) {
		plhs[0] = mxCreateStructMatrix(1, 1, knames, names);
		for( k=0; k<knames; k++ ) {
			mxSetFieldByNumber(plhs[0], 0, k, mxarrays[k]);;
		}
		mxFree(mxarrays);
	}
	if( knames ) {
		for( k=0; k<knames; k++ ) {
			mxFree(names[k]);
		}
		mxFree(names);
	}
	if( verbose ) {
		mexPrintf("\n");
	}
}

/**************************************************************************************************/

void fixvarname( char *name, char *newname, int verbose )
{
	char *oldn = name;
	char *newn = newname;
    unsigned char u;
	int badname = 0;

	if( !isalpha(*name) ) { /* Check that first char is a letter */
		badname = 1;
		*newname++ = 'A';
	}
	while( *name ) {
		if( isalnum(*name) || *name == '_' ) { /* Check for alphanumeric or underscore */
			*newname++ = *name++;
		} else {
			badname = 1;
            u = *name;
            if( u > 127 ) {
 			    u = upperASCII[u-128]; /* Replace bad characters with upper ASCII close matches */
            } else {
 			    u = '_'; /* Replace bad characters with underscores */
            }
	        *newname++ = u; /* Replace bad characters */
			name++;
		}
	}
	*newname = '\0'; /* Attach the final null termination character */
	if( verbose && badname ) {
		mexPrintf("Invalid Name:  <%s>\n",oldn);
		mexPrintf("Replaced With: <%s>\n",newn);
	}
}

/**************************************************************************************************/

void fixfieldnames(mxArray *mx, int verbose)
{
	int badname;
	mwSize i, j, n, f;
	char *fieldname, *fname, *fi, *fj, *newchar;
    unsigned char u;
	int z, nameclash;

	if( mx == NULL ) {
		; /* Nothing to fix */
	} else if( mxIsStruct(mx) ) {
		n = mxGetNumberOfElements(mx);
		f = mxGetNumberOfFields(mx);
		z = 0;
		for( j=0; j<f; j++ ) {
			badname = 0;
			fname = fieldname = (char *) mxGetFieldNameByNumber(mx,j);
			if( *fieldname == '\0' ) {
				mexWarnMsgTxt("Unable to fix empty fieldname.");
				continue;
			}
			if( !isalpha(*fieldname) ) { /* Check that first char is a letter */
				if( verbose ) {
					mexPrintf("Invalid Field Name:  <%s>\n",fname);
				}
				badname = 1;
				if( isdigit(*fieldname) ) {
					*fieldname = alphabet[*fieldname - '0'];
				} else {
					*fieldname = alphabet[z];
					z = (z+1) % 52;
				}
			}
			while( *fieldname ) {
				if( !isalnum(*fieldname) && *fieldname != '_' ) { /* Check for alphanumeric or underscore */
					if( !badname && verbose ) {
						mexPrintf("Invalid Field Name:  <%s>\n",fname);
					}
					badname = 1;
                    u = *fieldname;
                    if( u > 127 ) {
                        u = upperASCII[u-128]; /* Replace bad characters with close upper ASCII match */
                    } else {
					    u = '_'; /* Replace bad characters with underscores */
                    }
                    *fieldname++ = u; /* Replace bad characters */
				} else {
					fieldname++;
				}
			}
			if( verbose && badname ) {
				mexPrintf("Replaced With:       <%s>\n",fname);
			}
			nameclash = 1; /* Check for field name clashes */
			while( nameclash ) {
				nameclash = 0;
				for( i=0; i<f; i++ ) {
					fi = (char *) mxGetFieldNameByNumber(mx,i);
					for( j=i+1; j<f; j++ ) {
						fj = (char *) mxGetFieldNameByNumber(mx,j);
						if( fi && fj && strcmp(fi,fj) == 0 ) { /* Name clash */
							nameclash = 1;
							while( *fj ) { /* Fix the name clash */
							    newchar = alphanum;
								while( *newchar ) {
									if( *fj == *newchar ) {
										*fj = *(newchar+1);
										break;
									}
									newchar++;
								}
								fj++;
							}
						}
					}
				}
			}
		}
		for( i=0; i<n; i++ ) {
			for( j=0; j<f; j++ ) {
                fixfieldnames(mxGetFieldByNumber(mx,i,j), verbose);
			}
		}
	} else if( mxIsCell(mx) ) {
		n = mxGetNumberOfElements(mx);
		for( i=0; i<n; i++ ) {
            fixfieldnames(mxGetCell(mx,i), verbose);
		}
	}
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

/**************************************************************************************************/

int isamatch(const char *name, mxArray *compare)
{
	mwSize i, n;
	char *cname;

	if( compare == NULL ) return 1;
	if( mxIsChar(compare) ) {
		cname = mxArrayToString(compare);
//		i = strcmp(name,cname);
		i = !wildcmp(cname,name);
		mxFree(cname);
		return (i == 0);
	} else if( mxIsCell(compare) ) {
		n = mxGetNumberOfElements(compare);
		for( i=0; i<n; i++ ) {
			if( isamatch(name,mxGetCell(compare,i)) ) return 1;
		}
		return 0;
	} else {
		mexErrMsgTxt("Invalid cell. Must be char string with variable name.");
	}
}

/* wildcmp obtained from Jack Handy. Seems to work OK. */

int wildcmp(const char *wild, const char *str)
{
  const char *cp = NULL, *mp = NULL;

  while ((*str) && (*wild != '*')) {
    if ((*wild != *str) && (*wild != '?')) {
      return 0;
    }
    wild++;
    str++;
  }

  while (*str) {
    if (*wild == '*') {
      if (!*++wild) {
        return 1;
      }
      mp = wild;
      cp = str+1;
    } else if ((*wild == *str) || (*wild == '?')) {
      wild++;
      str++;
    } else {
      wild = mp;
      str = cp++;
    }
  }

  while (*wild == '*') {
    wild++;
  }
  return !*wild;
}
