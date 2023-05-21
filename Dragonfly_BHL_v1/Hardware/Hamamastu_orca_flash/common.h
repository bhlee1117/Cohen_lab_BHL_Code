/* * * * * * * * * * * * * * * * * * * *\
 *                                     *
 *   2018 Vicente Parot                *
 *   Cohen Lab - Harvard University    *
 *                                     *
 *   Modified from Hamamatsu DCAMAPI   *
 *                                     *
\* * * * * * * * * * * * * * * * * * * */

// modified from console/misc/common.h
//

#include "dcamapi4.h"
#include "dcamprop.h"

// #ifndef failed
// #define failed(error) (((int)error)<0) 
// #endif

#define MAX_MSG_LEN (260*8)//MAX_PATH
#define MAX_COMMAND 64
// #define round(x) (((x)>0)?(long)((x)+0.5):(long)((x)-0.5))

void dcamcon_show_dcamerr( HDCAM hdcam, DCAMERR errid, const char* apiname);

HDCAM dcamcon_init_open();
void dcamcon_show_dcamdev_info( HDCAM hdcam );
void message( const char* fmt );
void error( const char* fmt );
