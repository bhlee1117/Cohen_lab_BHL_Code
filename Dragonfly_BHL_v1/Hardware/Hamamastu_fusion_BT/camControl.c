/* * * * * * * * * * * * * * * * * * *\
 *                                    *
 *   2018 Vicente Parot               *
 *   Cohen Lab - Harvard University   *
 *                                    *
 \* * * * * * * * * * * * * * * * * * */

#include "camControl.h"

int sensor_h = 2304;
int sensor_w = 2304;

bool verboseFlag = true;
bool isCapturing = false;
bool isRecording = false;
uint16_t *lastImageData;
int lastImageDataWidth;
int lastImageDataHeight;
// DCAMBUF_FRAME *lastBufferFrame = NULL;
// DCAMBUF_FRAME *lastValidBufferFrame = NULL;

static DCAMERR 		err;
static HDCAM 		hdcam 	= NULL;
static HDCAMWAIT	hwait 	= NULL;
static HDCAMREC 	hrec 	= NULL;
static bool isDevOpen = false;
static bool isWaitAndBufferOpen = false;
static double exposureTime = .1;
static int32 hPos, hSize, vPos, vSize;
static int32 exposureTimeMilliseconds, readoutTimeMilliseconds;
double exposureTimeSeconds, readoutTimeSeconds;

char fpath_dcimg_copy[MAX_MSG_LEN];
char fext_dcimg_copy[MAX_MSG_LEN];
char fpath_bin_temp_copy[MAX_MSG_LEN];
char fext_bin_temp_copy[MAX_MSG_LEN];
char fpath_bin_tgt_copy[MAX_MSG_LEN];
char fext_bin_tgt_copy[MAX_MSG_LEN];

char msg[MAX_MSG_LEN];

extern HANDLE hConvertThread;
extern HANDLE ghMutexCapturing; 

void show_recording_status(void){
	// get recording status
	DCAMREC_STATUS recstatus;
	memset( &recstatus, 0, sizeof(recstatus) );
	recstatus.size	= sizeof(recstatus);
	
	err = dcamrec_status( hrec, &recstatus );
	if( failed(err) ){
		dcamcon_show_dcamerr( hdcam, err, "dcamrec_status()" );
	}else{
		sprintf(msg, "flags: 0x%08x, latest index: %06d, miss: %06d, total: %d\n", recstatus.flags, recstatus.currentframe_index, recstatus.missingframe_count, recstatus.totalframecount );
		message(msg);
	}
}

void dc_startup(void) {
	if (isDevOpen) {
		message("device already open\n");
	}else{
		// initialize DCAM-API and open device
		hdcam = dcamcon_init_open();
		if(!hdcam){ // failed open DCAM handle
			message("dcamcon_init_open returned NULL pointer\n");
			return;
		}
		message("dcamcon_init_open ok\n");
		isDevOpen = true;
		
		// show device information
		dcamcon_show_dcamdev_info( hdcam );
	}
	return;
}

bool dc_wait_buffer_open(int32 bufferFrames){
	if (isDevOpen) {
		// open wait handle
		DCAMWAIT_OPEN	waitopen;
		memset( &waitopen, 0, sizeof(waitopen) );
		waitopen.size	= sizeof(waitopen);
		waitopen.hdcam	= hdcam;

		err = dcamwait_open( &waitopen );
		if( failed(err) ) {
			dcamcon_show_dcamerr( hdcam, err, "dcamwait_open()" );
		}else{
			message("wait open\n");
			hwait = waitopen.hwait;

			// allocate buffer
			err = dcambuf_alloc( hdcam, bufferFrames );
			if( failed(err) ){
				dcamcon_show_dcamerr( hdcam, err, "dcambuf_alloc()" );
			}else{
				message("buffer allocated\n");
				return true;
			}
		}
	}else{
		message("no camera open\n");
	}
	return false;
}

bool dc_record_open_attach(int32 recordFrames, const char *fpath_dcimg, const char *fext_dcimg){
	if (isWaitAndBufferOpen) {
		strcpy(fpath_dcimg_copy,fpath_dcimg);
		strcpy(fext_dcimg_copy,fext_dcimg);
		
		// create file
		DCAMREC_OPEN	recopen;
		memset( &recopen, 0, sizeof(recopen) );
		recopen.size	= sizeof(recopen);
		recopen.path	= fpath_dcimg_copy;		// it should set new file name.
		recopen.ext		= fext_dcimg_copy;
		recopen.maxframepersession	= recordFrames;

		err = dcamrec_open( &recopen );
		if( failed(err) ){
			dcamcon_show_dcamerr( hdcam, err, "dcamrec_open()" );
		}else{
			message("record open\n");
			hrec = recopen.hrec;

			// attach recording handle to DCAM handle
			err = dcamcap_record( hdcam, hrec );
			if( failed(err) ){
				dcamcon_show_dcamerr( hdcam, err, "dcamcap_record()" );
			}else{
				message("record attached\n");
				return true;
			}
		}
	}else{
		message("buffer not allocated\n");
	}
	return false;
}

double get_property(const char idtxt[], int32 IDPROP){
	if (isDevOpen) {
		sprintf(msg,"get property %s ", idtxt);
		message(msg);
		double value;
		err = dcamprop_getvalue( hdcam, IDPROP, &value );
		if( failed(err) ){
			message("\n");
			dcamcon_show_dcamerr( hdcam, err, "dcamprop_getvalue()");
			return -1;
		}else{
			sprintf(msg, "returned value %g\n", value);
			message(msg);
			return value;
		}
	}else{
		message("no camera open\n");
	}
	return -1;
}

bool set_property(const char idtxt[], int32 IDPROP, double value){
	if (isDevOpen) {
		sprintf(msg,"set property %s with value %g\n", idtxt, value );
		message(msg);
		err = dcamprop_setvalue( hdcam, IDPROP, value );
		if( failed(err) )
		{
			dcamcon_show_dcamerr( hdcam, err, "dcamprop_setvalue()");
			get_property(idtxt, IDPROP);
			return false;
		}
	}else{
		message("no camera open\n");
	}
	return true;
}

bool set_subarray( int32 hpos, int32 hsize, int32 vpos, int32 vsize ){
	if(!isDevOpen){
		message("no camera open\n");
		return false;
	}
	set_property("SUBARRAYMODE, OFF", 	DCAM_IDPROP_SUBARRAYMODE, 	DCAMPROP_MODE__OFF);
	set_property("SUBARRAYHPOS", 		DCAM_IDPROP_SUBARRAYHPOS, 	hpos);
	set_property("SUBARRAYHSIZE", 		DCAM_IDPROP_SUBARRAYHSIZE, 	hsize);
	set_property("SUBARRAYVPOS", 		DCAM_IDPROP_SUBARRAYVPOS, 	vpos);
	set_property("SUBARRAYVSIZE", 		DCAM_IDPROP_SUBARRAYVSIZE, 	vsize);
	set_property("SUBARRAYMODE, ON", 	DCAM_IDPROP_SUBARRAYMODE, 	DCAMPROP_MODE__ON);
	return true;
}

bool set_centered_roi( int32 hSizeTry, int32 vSizeTry ){
	if(!isDevOpen){
		message("no camera open\n");
		return false;
	}
	hSize = round(hSizeTry/8.0)*8;
	vSize = round(vSizeTry/8.0)*8;
	if (hSize != hSizeTry){
		sprintf(msg,"rounding hSize from %d to %d\n", hSizeTry, hSize );
		message(msg);
	}
	if (vSize != vSizeTry){
		sprintf(msg,"rounding vSize from %d to %d\n", vSizeTry, vSize );
		message(msg);
	}
	hPos = (sensor_w-hSize)/2;
	vPos = (sensor_h-vSize)/2;
	set_subarray( hPos, hSize, vPos, vSize );
	return true;
}

bool set_arbitrary_roi( int32 hPosTry, int32 hSizeTry, int32 vPosTry, int32 vSizeTry ){
	if(!isDevOpen){
		message("no camera open\n");
		return false;
	}
	hSize = round(hSizeTry/4.0)*4;
	hSize = fmax(hSize,4);
	hSize = fmin(hSize,sensor_w);
	vSize = round(vSizeTry/4.0)*4;
	vSize = fmax(vSize,4);
	vSize = fmin(vSize,sensor_h);
	if (hSize != hSizeTry){
		sprintf(msg,"rounding hSize from %d to %d\n", hSizeTry, hSize );
		message(msg);
	}
	if (vSize != vSizeTry){
		sprintf(msg,"rounding vSize from %d to %d\n", vSizeTry, vSize );
		message(msg);
	}
	hPos = round(hPosTry/4.0)*4;
	hPos = fmax(hPos,0);
	hPos = fmin(hPos,sensor_w-hSize);
	vPos = round(vPosTry/4.0)*4;
	vPos = fmax(vPos,0);
	vPos = fmin(vPos,sensor_h-vSize);
	if (hPos != hPosTry){
		sprintf(msg,"rounding hPos from %d to %d\n", hPosTry, hPos );
		message(msg);
	}
	if (vPos != vPosTry){
		sprintf(msg,"rounding vPos from %d to %d\n", vPosTry, vPos );
		message(msg);
	}
	set_subarray( hPos, hSize, vPos, vSize );
	return true;
}

bool set_hsynctrigger( ) {
	if(!isDevOpen){
		message("no camera open\n");
		return false;
	}
	set_property("OUTPUTTRIGGER_KIND, PROGRAMABLE",	DCAM_IDPROP_OUTPUTTRIGGER_KIND, 	DCAMPROP_OUTPUTTRIGGER_KIND__PROGRAMABLE);
	set_property("OUTPUTTRIGGER_POLARITY, POSITIVE",DCAM_IDPROP_OUTPUTTRIGGER_POLARITY, DCAMPROP_OUTPUTTRIGGER_POLARITY__POSITIVE);
	set_property("OUTPUTTRIGGER_SOURCE, HSYNC", 	DCAM_IDPROP_OUTPUTTRIGGER_SOURCE, 	DCAMPROP_OUTPUTTRIGGER_SOURCE__HSYNC);
	set_property("OUTPUTTRIGGER_DELAY", 			DCAM_IDPROP_OUTPUTTRIGGER_DELAY, 	0);
	set_property("OUTPUTTRIGGER_PERIOD", 			DCAM_IDPROP_OUTPUTTRIGGER_PERIOD, 	2e-6);
	return true;
}

void aq_sync_stop() {
	if (isCapturing) {
		WaitForSingleObject(ghMutexCapturing, INFINITE);
			isCapturing = false;
		ReleaseMutex(ghMutexCapturing);
		dcamcap_stop( hdcam );
		message("capture stopped\n");
		// lastBufferFrame = NULL;
	}else{
		message("capture already stopped\n");
	}
	if (isRecording) {
		show_recording_status();
		dcamrec_close( hrec );

				// sprintf(msg,"converting file in background\n%s.%s\n",fpath_dcimg_copy,fext_dcimg_copy);
				// message(msg);
				// sprintf(msg,"to\n%s.%s\n",fpath_bin_temp_copy,fext_bin_temp_copy);
				// message(msg);

				unsigned int threadID;
				char msg[MAX_MSG_LEN];
				hConvertThread = (HANDLE)_beginthreadex( NULL, 0, &convertThreadFcn, NULL, 0, &threadID );
				message("converting file in background\n");
				sprintf(msg,"converting file in background\n%s.%s\n",fpath_dcimg_copy,fext_dcimg_copy);
				message(msg);


		message("record closed\n");
		isRecording = false;
	}else{
		message("record already stopped\n");
	}
	if (isWaitAndBufferOpen) {		
		dcambuf_release( hdcam );
		message("buffer released\n");
		dcamwait_close( hwait );
		message("wait closed\n");
		isWaitAndBufferOpen = false;
	}else{
		message("wait and buffer already closed\n");
	}
}

void aq_live_stop() {
	if (isCapturing) {
		WaitForSingleObject(ghMutexCapturing, INFINITE);
			isCapturing = false;
		ReleaseMutex(ghMutexCapturing);
		dcamcap_stop( hdcam );
		message("capture stopped\n");
		// lastBufferFrame = NULL;
	}else{
		message("capture already stopped\n");
	}
	// show_recording_status();
	// dcamrec_close( hrec );
	// message("record closed\n");
	// isRecording = false;
	if (isWaitAndBufferOpen) {
		dcambuf_release( hdcam );
		message("buffer released\n");
		dcamwait_close( hwait );
		message("wait closed\n");
		isWaitAndBufferOpen = false;
	}else{
		message("wait and buffer already closed\n");
	}
}

bool aq_live_restart(bool isCentered, double *inputROI, double binning, double exposureTime) {
	if(isDevOpen){
		if (isCapturing){
			// need to stop capturing
			if (isRecording){
				// synchronized recording to disk
				aq_sync_stop();
			}else{
				// live imaging
				aq_live_stop();
			} // both cases will have set isWaitAndBufferOpen to false also
		}else{
			// no need to stop capturing
		}
		// here assume capturing is stopped
		// set live imaging parameters
		set_property("BINNING",						DCAM_IDPROP_BINNING, 			binning);
		if (isCentered) {
			set_centered_roi(inputROI[0],inputROI[1]);
		} else {
			set_arbitrary_roi(inputROI[0],inputROI[1],inputROI[2],inputROI[3]);
		}
		set_property("READOUTSPEED, FASTEST",		DCAM_IDPROP_READOUTSPEED,		DCAMPROP_READOUTSPEED__FASTEST);
		set_property("TRIGGER_MODE, NORMAL",		DCAM_IDPROP_TRIGGER_MODE,		DCAMPROP_TRIGGER_MODE__NORMAL);
		set_property("TRIGGERSOURCE, INTERNAL",		DCAM_IDPROP_TRIGGERSOURCE,		DCAMPROP_TRIGGERSOURCE__INTERNAL);
		set_property("EXPOSURETIME",				DCAM_IDPROP_EXPOSURETIME,		exposureTime);
		double result;
		result = get_property("EXPOSURETIME",		DCAM_IDPROP_EXPOSURETIME);
		exposureTimeSeconds = result;
		exposureTimeMilliseconds = exposureTimeSeconds*1e3;
		result = get_property("TIMING_READOUTTIME",	DCAM_IDPROP_TIMING_READOUTTIME);
		readoutTimeSeconds = result;
		readoutTimeMilliseconds = readoutTimeSeconds*1e3;

		// start aq: wait and buffer, no recording, start camera
		if(!isWaitAndBufferOpen){
			int32 bufferFrames = 50;
			isWaitAndBufferOpen = dc_wait_buffer_open(bufferFrames);
			// if(!isRecording){
				// isRecording = dc_record_open_attach(recordFrames,"r:\\testfile","dcimg");
				// show_recording_status();
				if(!isCapturing){
					err = dcamcap_start( hdcam, DCAMCAP_START_SEQUENCE );
					if( failed(err) ){
						dcamcon_show_dcamerr( hdcam, err, "dcamcap_start()" );
					}else{
						message("capture started\n");
						WaitForSingleObject(ghMutexCapturing, INFINITE);
							isCapturing = true;
						ReleaseMutex(ghMutexCapturing);
						return true;
					}
				}else{
					message("capture is already started\n");
					// show_recording_status();
				}
			// }else{
				// message("record is already open\n");
				// show_recording_status();
			// }
		}else{
			message("wait and buffer are already open\n");
		}
	}else{
		message("no camera open\n");
	}	
	return false;
}

// bool aq_sync_prepare(double *inputROI, double exposureTime) {
bool aq_sync_prepare(bool isCentered, double *inputROI, double binning, double exposureTime) {
	if(isDevOpen){
		if (isCapturing){
			// need to stop capturing
			if (isRecording){
				// synchronized recording to disk
				aq_sync_stop();
			}else{
				// live imaging
				aq_live_stop();
			} // both cases will have set isWaitAndBufferOpen to false also
		}else{
			// no need to stop capturing
		}
		// here assume capturing is stopped
		// set sync imaging parameters
		set_property("BINNING",						DCAM_IDPROP_BINNING, 			binning);
		if (isCentered) {
			set_centered_roi(inputROI[0],inputROI[1]);
		} else {
			set_arbitrary_roi(inputROI[0],inputROI[1],inputROI[2],inputROI[3]);
		}
		set_property("READOUTSPEED, FASTEST",		DCAM_IDPROP_READOUTSPEED,		DCAMPROP_READOUTSPEED__FASTEST);
		set_property("TRIGGER_MODE, NORMAL",		DCAM_IDPROP_TRIGGER_MODE,		DCAMPROP_TRIGGER_MODE__NORMAL);
		set_property("TRIGGERSOURCE, EXTERNAL",		DCAM_IDPROP_TRIGGERSOURCE,		DCAMPROP_TRIGGERSOURCE__EXTERNAL);
		set_property("TRIGGERACTIVE, SYNCREADOUT",	DCAM_IDPROP_TRIGGERACTIVE,		DCAMPROP_TRIGGERACTIVE__SYNCREADOUT);
		set_property("TRIGGER_CONNECTOR, BNC",		DCAM_IDPROP_TRIGGER_CONNECTOR,	DCAMPROP_TRIGGER_CONNECTOR__BNC);
		set_property("TRIGGERPOLARITY, POSITIVE",	DCAM_IDPROP_TRIGGERPOLARITY,	DCAMPROP_TRIGGERPOLARITY__POSITIVE);
		//---- For Hamamatsu cameras with HSync firmware update only --------
        set_hsynctrigger();
        //---- For Hamamatsu cameras with HSync firmware update only --------
        
		set_property("EXPOSURETIME",				DCAM_IDPROP_EXPOSURETIME,		exposureTime);
		double result;
		result = get_property("EXPOSURETIME",		DCAM_IDPROP_EXPOSURETIME);
		exposureTimeSeconds = result;
		exposureTimeMilliseconds = exposureTimeSeconds*1e3;
		result = get_property("TIMING_READOUTTIME",	DCAM_IDPROP_TIMING_READOUTTIME);
		readoutTimeSeconds = result;
		readoutTimeMilliseconds = readoutTimeSeconds*1e3;
	}else{
		message("no camera open\n");
	}
	return false;	
}

bool aq_sync_start(int32 recordFrames, const char *fpath, const char *fext) {
	if(isDevOpen){
		if (isCapturing){
			// need to stop capturing
			if (isRecording){
				// synchronized recording to disk
				aq_sync_stop();
				
				// unsigned int threadID;
				// char msg[MAX_MSG_LEN];
				// hConvertThread = (HANDLE)_beginthreadex( NULL, 0, &convertThreadFcn, NULL, 0, &threadID );
				// message("converting file in background\n");
				// sprintf(msg,"converting file in background\n%s.%s\n",fpath_dcimg_copy,fext_dcimg_copy);
				// message(msg);
			}else{
				// live imaging
				aq_live_stop();
			} // both cases will have set isWaitAndBufferOpen to false also
		}else{
			// no need to stop capturing
		}
		if(!isWaitAndBufferOpen){
			int32 bufferFrames = 50;
			isWaitAndBufferOpen = dc_wait_buffer_open(bufferFrames);
			if(!isRecording){
				isRecording = dc_record_open_attach(recordFrames,fpath,fext);
				show_recording_status();
				if(!isCapturing){
					err = dcamcap_start( hdcam, DCAMCAP_START_SEQUENCE );
					if( failed(err) ){
						dcamcon_show_dcamerr( hdcam, err, "dcamcap_start()" );
					}else{
						message("capture started\n");
						WaitForSingleObject(ghMutexCapturing, INFINITE);
							isCapturing = true;
						ReleaseMutex(ghMutexCapturing);
						return true;
					}
				}else{
					message("capture is already started\n");
					show_recording_status();
				}
			}else{
				message("record is already open\n");
				show_recording_status();
			}
		}else{
			message("wait and buffer are already open\n");
		}
	}else{
		message("no camera open\n");
	}
	return false;
}

bool dc_shutdown(void) {
	if(isDevOpen){
		if (isCapturing){
			// need to stop capturing
			if (isRecording){
				// synchronized recording to disk
				aq_sync_stop();
				
				// unsigned int threadID;
				// char msg[MAX_MSG_LEN];
				// hConvertThread = (HANDLE)_beginthreadex( NULL, 0, &convertThreadFcn, NULL, 0, &threadID );
				// message("converting file in background\n");
				// sprintf(msg,"converting file in background\n%s.%s\n",fpath_dcimg_copy,fext_dcimg_copy);
				// message(msg);
			}else{
				// live imaging
				aq_live_stop();
			} // both cases will have set isWaitAndBufferOpen to false also
		}else{
			// no need to stop capturing
		}
		dcamdev_close( hdcam );
		message("dcamdev_close ok\n");
		// finalize DCAM-API
		dcamapi_uninit();
		message("dcamapi_uninit ok\n");
		isDevOpen = false;
		return true;
	}else{
		message("no camera open\n");
	}
	return false;
}


// double calc_histo( const void* buf, int32 rowbytes, DCAM_PIXELTYPE type, int32 width, int32 height )
// {
	// if( type != DCAM_PIXELTYPE_MONO16 )
	// {
		// // not implement
		// return -1;
	// }

	// // int32	cx = width / 10;
	// // int32	cy = height / 10;
	// // if( cx < 10 )	cx = 10;
	// // if( cy < 10 )	cy = 10;
	// // if( cx > width || cy > height )
	// // {
		// // // frame is too small
		// // return -1;
	// // }

	// // int32	ox = (width - cx) / 2;
	// // int32	oy = (width - cy) / 2;
	
	// int32	ox = 0;
	// int32	oy = 0;

	// const char*	src = (const char*)buf + rowbytes * oy;
	// double total = 0;

	// // calculate center sum
	// int32 x, y;
	// for( y=0; y < height; y++ )
	// {
		// const unsigned short*	s = (const unsigned short*)src + ox;
		// for( x = 0; x < width; x++ )
		// {
			// total += *s++;
		// }
	// }

	// return total / cx / cy;
// }

DCAMBUF_FRAME *aq_snap(){
	if (isCapturing){
		int32 nFrame = 1;

		// wait start param
		static DCAMWAIT_START	waitstart;
		memset( &waitstart, 0, sizeof(waitstart) );
		waitstart.size		= sizeof(waitstart);
		waitstart.eventmask	= DCAMWAIT_CAPEVENT_FRAMEREADY;
		waitstart.timeout	= exposureTimeMilliseconds + readoutTimeMilliseconds; // milliseconds

		// prepare frame param
		static DCAMBUF_FRAME	bufframe;
		memset( &bufframe, 0, sizeof(bufframe) );
		bufframe.size		= sizeof(bufframe);
		bufframe.iFrame		= -1;				// latest captured image

		DCAMERR err;

		// wait image
		err = dcamwait_start( hwait, &waitstart );
		if( failed(err) ){
			dcamcon_show_dcamerr( hdcam, err, "dcamwait_start()" );
			return NULL;
		}
		message("wait started\n");

		// access image
		err = dcambuf_lockframe( hdcam, &bufframe );
		if( failed(err) ){
			dcamcon_show_dcamerr( hdcam, err, "dcambuf_lockframe()" );
			return NULL;
		}
		message("frame locked\n");
		
		return &bufframe;
	}else{
		message("camera is not capturing\n");
	}
	return NULL;
}

DCAMBUF_FRAME *aq_thread_snap(){
	// if (isCapturing){ // camera is always capturing when this function is called
	int32 nFrame = 1;

	// wait start param
	DCAMWAIT_START	waitstart;
	memset( &waitstart, 0, sizeof(waitstart) );
	waitstart.size		= sizeof(waitstart);
	waitstart.eventmask	= DCAMWAIT_CAPEVENT_FRAMEREADY;
	waitstart.timeout	= 2*(exposureTimeMilliseconds + readoutTimeMilliseconds); // milliseconds

	// prepare frame param
	static DCAMBUF_FRAME	bufframe;
	memset( &bufframe, 0, sizeof(bufframe) );
	bufframe.size		= sizeof(bufframe);
	bufframe.iFrame		= -1;				// latest captured image

	DCAMERR err;

	// wait image
	err = dcamwait_start( hwait, &waitstart );
	if( failed(err) ){
		// dcamcon_show_dcamerr( hdcam, err, "dcamwait_start()" );
		return NULL;
	}

	// access image
	err = dcambuf_lockframe( hdcam, &bufframe );
	if( failed(err) ){
		// dcamcon_show_dcamerr( hdcam, err, "dcambuf_lockframe()" );
		return NULL;
	}
	return &bufframe;
}

// void dc_test() {
	// // set binning value to the camera
	// set_property("BINNING", DCAM_IDPROP_BINNING, DCAMPROP_BINNING__1);

// }

// void __mexFunction__camControl( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){

	// char *input_command;
	// int buflen;
	// double *centROI;
	// double expTime;


	// // message("called with %d input argument%s\n",nrhs,(nrhs-1)?"s":"");
    // if(nrhs<1){
		// error("for startup, use the startup command\n");
		// return;
	// }else{
		// if ( mxIsChar(prhs[0]) != 1){
			// error("first input must be a string: camControl('startup')\n");
			// return;
		// }

		// if (mxGetM(prhs[0])!=1){
			// error("first input must be a row vector string: camControl('startup')\n");
			// return;
		// }
		
		// // get the length of the input string
		// buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
		
		// // copy the string data from prhs[0] into a C string input_command
		// input_command = mxArrayToString(prhs[0]);
	// }		
    
    // if(input_command == NULL){
		// error("could not convert input to string\n");
	// }
	
	// if (!strcmp(input_command, "startup")) { 					////  startup
		// message("running startup\n");
		// dc_startup();
	// } else if (!strcmp(input_command, "shutdown")) { 			////  shutdown
		// message("running shutdown\n");
		// dc_shutdown();
	// } else if (!strcmp(input_command, "aq_sync_prepare")) { 	////  aq_sync_prepare
		// message("running aq_sync_prepare\n");
		// if( (nrhs!=3) || 
			// (!mxIsNumeric(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=2 || mxGetNumberOfDimensions(prhs[1])!=2) || // 1st parameter, ROI
			// (!mxIsNumeric(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1 || mxGetNumberOfDimensions(prhs[2])!=2) ) { // 2nd parameter, exposure time
			// error("aq_sync_prepare requires two parameters: camControl('aq_sync_prepare',inputROI,exposureTime)\n");
			// error("inputROI must be a 1x2 vector and exposureTime must be a numeric scalar\n");
			// return;
		// }
		// centROI = ((double *)(mxGetData(prhs[1])));
		// expTime = *((double *)(mxGetData(prhs[2])));
		// aq_sync_prepare(centROI,expTime);
	// } else if (!strcmp(input_command, "aq_sync_start")) { 		////  aq_sync_start
		// message("running aq_sync_start\n");
		// if( (nrhs!=2) ||
			// (!mxIsNumeric(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1 || mxGetNumberOfDimensions(prhs[1])!=2) ) {
			// error("aq_sync_start requires one parameters: camControl('aq_sync_start',recordFrames)\n");
			// error("recordFrames must be a numeric scalar\n");
			// return;
		// }
		// double recFrames = *((double *)(mxGetData(prhs[1])));
		// isCapturing = aq_sync_start(recFrames);
	// } else if (!strcmp(input_command, "test")) { 				////  test
		// message("testing nothing\n");
	// } else if (!strcmp(input_command, "aq_sync_stop")) { 		////  aq_sync_stop
		// message("running aq_sync_stop\n");
		// aq_sync_stop();
	// } else if (!strcmp(input_command, "aq_live_restart")) { 	////  aq_live_restart
		// message("running aq_live_restart\n");
		// if( (nrhs!=3) || 
			// (!mxIsNumeric(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=2 || mxGetNumberOfDimensions(prhs[1])!=2) || // 1st parameter, ROI
			// (!mxIsNumeric(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1 || mxGetNumberOfDimensions(prhs[2])!=2) ) { // 2nd parameter, exposure time
			// error("aq_live_restart requires two parameters: camControl('aq_live_restart',inputROI,exposureTime)\n");
			// error("inputROI must be a 1x2 vector and exposureTime must be a numeric scalar\n");
			// return;
		// }
		// centROI = ((double *)(mxGetData(prhs[1])));
		// expTime = *((double *)(mxGetData(prhs[2])));
		// isCapturing = aq_live_restart(centROI,expTime);
	// } else if (!strcmp(input_command, "aq_snap")) { 			////  aq_snap
		// message("running aq_snap\n");
		// if (nrhs>1){
			// error("aq_snap requires at most one output: img = camControl('aq_snap')\n");
		// }
		// if (aq_snap()) {
			// // there is a new image
			// mxArray *mxOut = mxCreateUninitNumericMatrix(lastBufferFrame->width, lastBufferFrame->height, mxUINT16_CLASS, mxREAL);
			// unsigned short *outBuf = (unsigned short *)mxGetData(mxOut);
			// memcpy(outBuf, lastBufferFrame->buf, lastBufferFrame->height*lastBufferFrame->width*sizeof(unsigned short));
			// sprintf(msg,"%d %d\n", lastBufferFrame->height, lastBufferFrame->width);
			// printf(msg);
			
			// plhs[0] = mxOut;
			// message("return new image\n");
		// }else{
			// if (lastBufferFrame) {
				// // the image is old but still a valid pointer
				// mxArray *mxOut = mxCreateUninitNumericMatrix(lastBufferFrame->width, lastBufferFrame->height, mxUINT16_CLASS, mxREAL);
				// unsigned short *outBuf = (unsigned short *)mxGetData(mxOut);
				// memcpy(outBuf, lastBufferFrame->buf, lastBufferFrame->height*lastBufferFrame->width*sizeof(unsigned short));
				// plhs[0] = mxOut;
				// message("return old image\n");
			// }else{
				// mxArray *mxOut = mxCreateUninitNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
				// unsigned short *outBuf = (unsigned short *)mxGetData(mxOut);
				// outBuf[0] = 0;
				// plhs[0] = mxOut;
				// message("return no image\n");
			// }
		// }
	// } else if (!strcmp(input_command, "aq_live_stop")) { 		////  aq_live_stop
		// message("running aq_live_stop\n");
		// aq_live_stop();
	// } else if (!strcmp(input_command, "verbose")) { 		////  aq_live_stop
		// verboseFlag = !verboseFlag;
		// sprintf(msg,"verbose flag set to %s\n", verboseFlag?"true":"false");
		// printf(msg);
	// } else { 													////  unknown
		// sprintf(msg,"unknown command '%s'\n", input_command);
		// message(msg);
	// }
	// return;
  // }