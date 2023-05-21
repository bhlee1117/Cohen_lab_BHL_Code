/* * * * * * * * * * * * * * * * * * *\
 *                                    *
 *   2018 Vicente Parot               *
 *   Cohen Lab - Harvard University   *
 *                                    *
 \* * * * * * * * * * * * * * * * * * */

#include "cameraMex.h"

extern bool verboseFlag;
extern bool isCapturing;
extern uint16_t *lastImageData;
extern int lastImageDataWidth;
extern int lastImageDataHeight;
// extern DCAMBUF_FRAME *lastBufferFrame;
// extern DCAMBUF_FRAME *lastValidBufferFrame;
extern HANDLE ghMutexFrame; 
extern HANDLE ghMutexCapturing; 

extern char fpath_bin_temp_copy[MAX_MSG_LEN];
extern char fext_bin_temp_copy[MAX_MSG_LEN];
extern char fpath_bin_tgt_copy[MAX_MSG_LEN];
extern char fext_bin_tgt_copy[MAX_MSG_LEN];

extern double exposureTimeSeconds, readoutTimeSeconds;

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){

	char *input_command, *param1, *param2, *param3, *param4, *param5, *param6;
	int buflen;
	double *inputROI, binning, expTime;
	bool isROICentered;
	char msg[MAX_MSG_LEN];
	DCAMBUF_FRAME *newFrame = NULL;

	// message("called with %d input argument%s\n",nrhs,(nrhs-1)?"s":"");
    if(nrhs<1){
		error("for startup, use the startup command\n");
		return;
	}else{
		if ( mxIsChar(prhs[0]) != 1){
			error("first input must be a string: cameraMex('startup')\n");
			return;
		}

		if (mxGetM(prhs[0])!=1){
			error("first input must be a row vector string: cameraMex('startup')\n");
			return;
		}
		
		// get the length of the input string
		buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
		
		// copy the string data from prhs[0] into a C string input_command
		input_command = mxArrayToString(prhs[0]);
	}		
    
    if(input_command == NULL){
		error("could not convert input to string\n");
	}
	
	if (!strcmp(input_command, "startup")) { 					////  startup
		uint16_t *windowParams = NULL;
		switch (nrhs) {
			case 1:
				break;
			case 6:
				windowParams = (uint16_t*)malloc(5*sizeof(uint16_t));
				windowParams[0] = (uint16_t)*((double *)(mxGetData(prhs[1])));
				windowParams[1] = (uint16_t)*((double *)(mxGetData(prhs[2])));
				windowParams[2] = (uint16_t)*((double *)(mxGetData(prhs[3])));
				windowParams[3] = (uint16_t)*((double *)(mxGetData(prhs[4])));
				windowParams[4] = (uint16_t)*((double *)(mxGetData(prhs[5])));
				break;
			default:
				error("startup must have 0 or 5 additional parameters\n");
				return;
				break;
		}
		message("running startup\n");
		load_test_image();
		dc_startup(); // dcam startup
		if (mexIsLocked()){
			mexErrMsgTxt("the program is already open, close and open again");
		}
		start_SDL_window(windowParams);
		mexLock(); // lock mex so that no-one accidentally clears mex
		launch_threads();
	} else if (!strcmp(input_command, "shutdown")) { 			////  shutdown
		message("running shutdown\n");
		dc_shutdown();

		// make sure that the thread was started using "mex locked" status
		if (!mexIsLocked()){
			mexErrMsgTxt("Thread not initialized yet."); // This function will return control to MATLAB
		}
		close_SDL_window_and_threads();

		mexUnlock();
		mexPrintf( "mex unlocked ok\n" );
	} else if (!strcmp(input_command, "aq_sync_prepare")) { 	////  aq_sync_prepare
		message("running aq_sync_prepare\n");
		if( (nrhs!=4) || 
			(!mxIsNumeric(prhs[1]) || !(mxGetNumberOfElements(prhs[1])==2 || mxGetNumberOfElements(prhs[1])==4) || mxGetNumberOfDimensions(prhs[1])!=2) || // 1st parameter, ROI
			(!mxIsNumeric(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1 || mxGetNumberOfDimensions(prhs[2])!=2) || // 2nd parameter, binning
			(!mxIsNumeric(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1 || mxGetNumberOfDimensions(prhs[3])!=2) ) { // 3nd parameter, exposure time
			error("aq_sync_prepare requires four parameters: cameraMex('aq_live_restart',centeredROI,binning,exposureTime)\n");
			error("centeredROI must be a 2 or 4-element vector\n");
			error("binning must be a numeric scalar with value 1, 2, 4\n");
			error("exposureTime must be a numeric scalar\n");
			return;
		}
		inputROI = ((double *)(mxGetData(prhs[1])));
		binning = *((double *)(mxGetData(prhs[2])));
		expTime = *((double *)(mxGetData(prhs[3])));
		isROICentered = !(2-mxGetNumberOfElements(prhs[1]));
		switch ((int)binning) {
			case 1:
				aq_sync_prepare(isROICentered, inputROI, DCAMPROP_BINNING__1, expTime);
				break;
			case 2:
				aq_sync_prepare(isROICentered, inputROI, DCAMPROP_BINNING__2, expTime);
				break;
			case 4:
				aq_sync_prepare(isROICentered, inputROI, DCAMPROP_BINNING__4, expTime);
				break;
			default:
				error("binning must be a numeric scalar with value 1, 2, or 4\n");
				return;
				break;
		}
		plhs[0] = mxCreateDoubleScalar(exposureTimeSeconds);
		plhs[1] = mxCreateDoubleScalar(readoutTimeSeconds);
// 		if( (nrhs!=3) || 
// 			(!mxIsNumeric(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=2 || mxGetNumberOfDimensions(prhs[1])!=2) || // 1st parameter, ROI
// 			(!mxIsNumeric(prhs[2]) || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1 || mxGetNumberOfDimensions(prhs[2])!=2) ) { // 2nd parameter, exposure time
// 			error("aq_sync_prepare requires two parameters: cameraMex('aq_sync_prepare',centeredROI,exposureTime)\n");
// 			error("centeredROI must be a 1x2 vector and exposureTime must be a numeric scalar\n");
// 			return;
// 		}
// 		inputROI = ((double *)(mxGetData(prhs[1])));
// 		expTime = *((double *)(mxGetData(prhs[2])));
// 		aq_sync_prepare(inputROI,expTime);
	} else if (!strcmp(input_command, "aq_sync_start")) { 		////  aq_sync_start
		message("running aq_sync_start\n");
		if( (nrhs!=8) ||
			(!mxIsNumeric(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1 || mxGetNumberOfDimensions(prhs[1])!=2) ) {
			error("aq_sync_start requires three parameters: cameraMex('aq_sync_start',recordFrames,\n");
			error("    dcimg_fpath,dcimg_fext,\n");
			error("    bin_temp_fpath,bin_temp_fext,\n");
			error("    bin_tgt_fpath,bin_tgt_fext)\n");
			error("recordFrames must be a numeric scalar\n");
			error("fpath must be double back-slashed path strings\n");
			error("fext must be extension strings\n");
			return;
		}
		double recFrames = *((double *)(mxGetData(prhs[1])));
		param1 = mxArrayToString(prhs[2]);
		param2 = mxArrayToString(prhs[3]);
		param3 = mxArrayToString(prhs[4]);
		param4 = mxArrayToString(prhs[5]);
		param5 = mxArrayToString(prhs[6]);
		param6 = mxArrayToString(prhs[7]);
		strcpy(fpath_bin_temp_copy,param3);
		strcpy(fext_bin_temp_copy,param4);
		strcpy(fpath_bin_tgt_copy,param5);
		strcpy(fext_bin_tgt_copy,param6);
		aq_sync_start(recFrames,param1,param2);
		message("aq_sync_start ok\n");
	} else if (!strcmp(input_command, "test")) { 				////  test
		message("testing nothing\n");
	} else if (!strcmp(input_command, "aq_sync_stop")) { 		////  aq_sync_stop
		message("running aq_sync_stop\n");
		aq_sync_stop();
	} else if (!strcmp(input_command, "aq_live_restart")) { 	////  aq_live_restart
		message("running aq_live_restart\n");
		if( (nrhs!=4) || 
			(!mxIsNumeric(prhs[1]) || !(mxGetNumberOfElements(prhs[1])==2 || mxGetNumberOfElements(prhs[1])==4) || mxGetNumberOfDimensions(prhs[1])!=2) || // 1st parameter, ROI
			(!mxIsNumeric(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1 || mxGetNumberOfDimensions(prhs[2])!=2) || // 2nd parameter, binning
			(!mxIsNumeric(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1 || mxGetNumberOfDimensions(prhs[3])!=2) ) { // 3nd parameter, exposure time
			error("aq_live_restart requires two parameters: cameraMex('aq_live_restart',centeredROI,binning,exposureTime)\n");
			error("centeredROI must be a 2 or 4-element vector\n");
			error("binning must be a numeric scalar with value 1, 2, 4\n");
			error("exposureTime must be a numeric scalar\n");
			return;
		}
		inputROI = ((double *)(mxGetData(prhs[1])));
		binning = *((double *)(mxGetData(prhs[2])));
		expTime = *((double *)(mxGetData(prhs[3])));
		// aq_live_restart(inputROI,expTime);
		isROICentered = !(2-mxGetNumberOfElements(prhs[1]));
		switch ((int)binning) {
			case 1:
				aq_live_restart(isROICentered, inputROI, DCAMPROP_BINNING__1, expTime);
				break;
			case 2:
				aq_live_restart(isROICentered, inputROI, DCAMPROP_BINNING__2, expTime);
				break;
			case 4:
				aq_live_restart(isROICentered, inputROI, DCAMPROP_BINNING__4, expTime);
				break;
			default:
				error("binning must be a numeric scalar with value 1, 2, or 4\n");
				return;
				break;
		}
		plhs[0] = mxCreateDoubleScalar(exposureTimeSeconds);
		plhs[1] = mxCreateDoubleScalar(readoutTimeSeconds);
	} else if (!strcmp(input_command, "aq_snap")) { 			////  aq_snap
		message("running aq_snap\n");
		if (nrhs>1){
			error("aq_snap requires at most one output: img = cameraMex('aq_snap')\n");
		}
		newFrame = aq_snap();
		if (newFrame) {
			WaitForSingleObject(ghMutexFrame, INFINITE);
				memcpy(lastImageData, newFrame->buf, newFrame->width*newFrame->height*sizeof(uint16_t));
				lastImageDataWidth  = newFrame->width;
				lastImageDataHeight = newFrame->height;
			ReleaseMutex(ghMutexFrame);
			message("got new image, returning\n");
		}else{
			error("timeout with no new image, returning last\n");
		}

		// there is a new image
		mxArray *mxOut = mxCreateUninitNumericMatrix(lastImageDataWidth, lastImageDataHeight, mxUINT16_CLASS, mxREAL);
		uint16_t *outBuf = (uint16_t *)mxGetData(mxOut);
		memcpy(outBuf, lastImageData, lastImageDataWidth*lastImageDataHeight*sizeof(uint16_t));
		// sprintf(msg,"dimensions: %d %d\n", lastImageDataWidth, lastImageDataHeight);
		// message(msg);
		plhs[0] = mxOut;
	// } else if (!strcmp(input_command, "aq_thread_snap")) { 			////  aq_thread_snap test
		// // message("running aq_snap\n");
		// if (nrhs>1){
			// error("aq_thread_snap requires at most one output: img = cameraMex('aq_thread_snap')\n");
		// }
		// message("start 10 snaps\n");
		// for (int num = 0; num < 10; num++){
			// newFrame = aq_thread_snap();
			// if (newFrame) {
				// message("new frame acquired\n");
				// WaitForSingleObject(ghMutexFrame, INFINITE);
					// memcpy(lastImageData, newFrame->buf, lastImageDataWidth*lastImageDataHeight*sizeof(uint16_t));
					// lastImageDataWidth  = newFrame->width;
					// lastImageDataHeight = newFrame->height;
				// ReleaseMutex(ghMutexFrame);
				// message("new frame copied\n");
			// }else{
				// message("null frame\n");
			// }
		// }
		// message("end 10 snaps\n");
	} else if (!strcmp(input_command, "aq_live_stop")) { 		////  aq_live_stop
		message("running aq_live_stop\n");
		aq_live_stop();
	} else if (!strcmp(input_command, "verbose")) { 		////  aq_live_stop
		verboseFlag = !verboseFlag;
		sprintf(msg,"verbose flag set to %s\n", verboseFlag?"true":"false");
		printf(msg);
	} else { 													////  unknown
		sprintf(msg,"unknown command '%s'\n", input_command);
		message(msg);
	}
	return;
  }