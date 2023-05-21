/* * * * * * * * * * * * * * * * * * * *\
 *                                     *
 *   2018 Vicente Parot                *
 *   Cohen Lab - Harvard University    *
 *                                     *
 *   Modified from Hamamatsu DCIMGAPI  *
 *                                     *
\* * * * * * * * * * * * * * * * * * * */

#ifndef _include_dcimg2bin_h_
#define _include_dcimg2bin_h_

#include "mex.h"
#include <windows.h>
#include "dcimgapi.h"

#include "common.h"

#ifndef mex_h
#define myErrMsgFcn(MSG) {printf(MSG);printf("\n");return;}
#else
#define myErrMsgFcn(MSG) {mexErrMsgTxt(MSG);}
#endif

// #ifndef failed
// #define failed(error) (((int)error)<0) 
// #endif

void dcimgcon_show_dcimgerr( DCIMG_ERR errid, const char* apiname);
HDCIMG dcimgcon_init_open( const char* filename );
char* get_extension( char* p );
bool convert_file(const char *dcimg_fname,const char *fname_bin_temp);

#endif // define _include_dcimg2bin_h_