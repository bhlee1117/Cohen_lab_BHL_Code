
#ifndef pvcamControl_h_
#define pvcamControl_h_


#include <stdio.h>
#include <cmath>
#include "Common.h"


typedef struct BUF_FRAME
{
	void* buf{nullptr};
	uint32_t width{0};
	uint32_t height{0};
} BUF_Frame;

BUF_FRAME *aq_snap(void);
bool pvcam_startup(void);
void pixel_filter_off(void);
void output_pp_parameter_state(std::string &out_msg);
void output_speedtable(std::string &out_msg);
double get_property(const char idtxt[], int32 IDPROP);
bool set_property(const char idtxt[], int32 IDPROP, double value);
bool set_subarray( int32 hpos, int32 hsize, int32 vpos, int32 vsize );
bool set_centered_roi( double hSizeTry, double vSizeTry );
bool set_arbitrary_roi( double hPosTry, double hSizeTry, double vPosTry, double vSizeTry );
bool aq_live_stop();
bool aq_sync_prepare(bool isCentered, double *inputROI, double binning, double exposureTime);
bool aq_sync_start_Fcn(int32 recordFrames, const char *fpath, const char *fext);
void aq_sync_start(int32 recordFrames, const char *fpath, const char *fext);
bool pvcam_shutdown(void);
bool aq_live_restart(bool isCentered, double *centeredROI, double binning, double exposureTime);
bool wait_for_thread_abort(CameraContext* ctx, uns32 timeoutMs);
bool stream_data_to_file(char* fpath, uns8* circBufferPtr, const uns32 circBufferSize, int totalFrames);


#endif // define _include_pvcamControl_h_
