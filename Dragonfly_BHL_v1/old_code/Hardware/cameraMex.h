/* * * * * * * * * * * * * * * * * * *\
 *                                    *
 *   2018 Vicente Parot               *
 *   Cohen Lab - Harvard University   *
 *                                    *
 \* * * * * * * * * * * * * * * * * * */

#ifndef _include_cameraMex_h_
#define _include_cameraMex_h_

#include <SDL.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"
#include <windows.h> // for threads, mutex
#include <process.h> // for threads, mutex
#include "dcamapi4.h"
#include "dcamprop.h"
#include "common.h"

#include "camControl.h"
#include "graphicThreads.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );

#endif // define _include_cameraMex_h_