/* * * * * * * * * * * * * * * * * * *\
 *                                    *
 *   2018 Vicente Parot               *
 *   Cohen Lab - Harvard University   *
 *                                    *
 \* * * * * * * * * * * * * * * * * * */

#ifndef _include_camControl_h_
#define _include_camControl_h_

#include <SDL.h>
#include <stdio.h>
#include <math.h>
// #include "mex.h"
// #include "matrix.h"
#include <windows.h> // for threads, mutex
#include <process.h> // for threads, mutex
#include "dcamapi4.h"
#include "dcamprop.h"
#include "common.h"

#include "dcimg2bin.h"
#include "graphicThreads.h"

void show_recording_status(void);
void dc_startup(void);
bool dc_wait_buffer_open(int32 bufferFrames);
bool dc_record_open_attach(int32 recordFrames, const char *fpath, const char *fext);
double get_property(const char idtxt[], int32 IDPROP);
bool set_property(const char idtxt[], int32 IDPROP, double value);
bool set_subarray( int32 hpos, int32 hsize, int32 vpos, int32 vsize );
bool set_centered_roi( int32 hSizeTry, int32 vSizeTry );
bool set_arbitrary_roi( int32 hPosTry, int32 hSizeTry, int32 vPosTry, int32 vSizeTry );
bool set_hsynctrigger( );
void aq_sync_stop();
void aq_live_stop();
bool aq_sync_prepare(bool isCentered, double *inputROI, double binning, double exposureTime);
bool aq_sync_start(int32 recordFrames, const char *fpath, const char *fext);
bool dc_shutdown(void);
bool aq_live_restart(bool isCentered, double *centeredROI, double binning, double exposureTime);
DCAMBUF_FRAME *aq_snap(void);
DCAMBUF_FRAME *aq_thread_snap(void);

#endif // define _include_camControl_h_
