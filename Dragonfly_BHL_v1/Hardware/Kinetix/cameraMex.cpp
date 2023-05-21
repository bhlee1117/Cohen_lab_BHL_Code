#include "cameraMex.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;
using namespace matlab::engine;

std::shared_ptr<MATLABEngine> matlabPtr;

extern bool verboseFlag;
extern bool isCapturing;
extern uint16_t lastImageData[IMG_WIDTH * IMG_HEIGHT];
extern uint16_t lastImageDataWidth;
extern uint16_t lastImageDataHeight;

extern std::mutex ghMutexFrame; 
extern std::mutex ghMutexCapturing; 

extern std::string out_pp_param;
extern std::string out_speedtable;

extern char fpath_dcimg_copy[MAX_MSG_LEN];
extern char fext_dcimg_copy[MAX_MSG_LEN];
extern char fpath_bin_temp_copy[MAX_MSG_LEN];
extern char fext_bin_temp_copy[MAX_MSG_LEN];
extern char fpath_bin_tgt_copy[MAX_MSG_LEN];
extern char fext_bin_tgt_copy[MAX_MSG_LEN];

extern double exposureTimeSeconds, readoutTimeSeconds;

//! Extracts the pointer to underlying data from the non-const iterator (`TypedIterator<T>`).
/*! This function does not throw any exceptions. */
template <typename T>
inline T* toPointer(const matlab::data::TypedIterator<T>& it) MW_NOEXCEPT {
    static_assert(std::is_arithmetic<T>::value && !std::is_const<T>::value,
        "Template argument T must be a std::is_arithmetic and non-const type.");
    return it.operator->();
}

template <typename T>
inline T* getPointer(matlab::data::TypedArray<T>& arr) MW_NOEXCEPT {
    static_assert(std::is_arithmetic<T>::value, "Template argument T must be a std::is_arithmetic type.");
    return toPointer(arr.begin());
}
template <typename T>
inline const T* getPointer(const matlab::data::TypedArray<T>& arr) MW_NOEXCEPT {
    return getPointer(const_cast<matlab::data::TypedArray<T>&>(arr));
}

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
		ArrayFactory factory;
		const char *param1, *param2, *param3, *param4, *param5, *param6;
		int buflen;
		double *inputROI, binning, expTime;
		bool isROICentered;
		BUF_FRAME *newFrame = NULL;
		
		matlabPtr = getEngine();
		
        std::string input_command = CharArray(inputs[0]).toAscii();
		
        if(inputs.size() < 1){
				matlabPtr->eval(u"error('for startup, use the startup command\\n')");
				return;
			}else{
				if ( inputs[0].getType() != ArrayType::CHAR){
					matlabPtr->eval(u"error('first input must be a string: cameraMex('startup')\\n')");
					return;
				}

				if (inputs[0].getDimensions()[0] != 1){
					matlabPtr->eval(u"error('first input must be a row vector string: cameraMex('startup')\\n')");
					return;
				}
				
				// get the length of the input string
				auto buflen = inputs[0].getNumberOfElements();
				
				// copy the string data from inputs[0] into a C++ string input_command
// 				*input_command = *(CharArray(inputs[0]).toAscii().c_str());
                
			}		
			
// 			if(input_command == NULL){
// 				matlabPtr->eval(u"error('could not convert input to string\n'");
// 			}
			
			if (!std::strcmp(input_command.c_str(), "startup")) { 					////  startup
				uint16_t *windowParams = NULL;
				switch (inputs.size()) {
					case 1:
						break;
					case 6:
                    {
// 						windowParams = (uint16_t*)malloc(5*sizeof(uint16_t));
                        windowParams = (uint16_t*)malloc(5*sizeof(double));
                        
//                         matlabPtr->eval(u"fprintf('input_command conversion \\n')");
                        const TypedArray<double> input1 = std::move(inputs[1]);
						windowParams[0] = input1[0];
                        
                        
                        const TypedArray<double> input2 = std::move(inputs[2]);
                        windowParams[1] = input2[0];
                        
                        const TypedArray<double> input3 = std::move(inputs[3]);
						windowParams[2] = input3[0];
                        
                        const TypedArray<double> input4 = std::move(inputs[4]);
						windowParams[3] = input4[0];
                        
                        const TypedArray<double> input5 = std::move(inputs[5]);
						windowParams[4] = input5[0];
                        
						break;
                    }
					default:
						matlabPtr->eval(u"error('startup must have 0 or 5 additional parameters\\n')");
						return;
						break;
				}
				matlabPtr->eval(u"fprintf('running startup\\n')");
				load_test_image();
				if(pvcam_startup()){ // pvcam startup
                    start_SDL_window(windowParams);
                    mexLock(); // lock mex so that no-one accidentally clears mex
                    launch_threads();
                    CharArray out = factory.createCharArray("");
//                     out = factory.createCharArray(out_pp_param);
//                     matlabPtr->feval(u"fprintf",out);
//                     out = factory.createCharArray(out_speedtable);
//                     matlabPtr->feval(u"fprintf",out);
//                     
//                     output_speedtable(out_msg);
//                     out = factory.createCharArray(out_msg);
//                     matlabPtr->feval(u"fprintf",out);
                    
                    matlabPtr->eval(u"fprintf('startup...ok\\n')");
                }
			} else if (!std::strcmp(input_command.c_str(), "shutdown")) { 			////  shutdown
				matlabPtr->eval(u"fprintf('running shutdown\\n')");
				close_SDL_window_and_threads();
                if(pvcam_shutdown()){
                    matlabPtr->eval(u"fprintf('PVCAM shutdown ok \\n')");
                }
				// make sure that the thread was started using "mex locked" status
                std::this_thread::sleep_for(std::chrono::milliseconds(500));
				mexUnlock();
				matlabPtr->eval(u"fprintf('mex unlocked ok \\n')");
			} else if (!std::strcmp(input_command.c_str(), "aq_sync_prepare")) { 	////  aq_sync_prepare
				matlabPtr->eval(u"fprintf('running aq_sync_prepare\\n')");
// 				if( (inputs.size() != 4) || 
// 					((!inputs[1].getNumberOfElements() == 2 || inputs[1].getNumberOfElements() == 4) || getNumElements(inputs[1].getDimensions()) != 2) || // 1st parameter, ROI
// 					(inputs[2].getNumberOfElements() != 1 || getNumElements(inputs[2].getDimensions()) != 2) || // 2nd parameter, binning
// 					(inputs[3].getNumberOfElements() !=1 || getNumElements(inputs[3].getDimensions()) != 2) ) { // 3nd parameter, exposure time
// 					matlabPtr->eval(u"error('aq_sync_prepare requires four parameters: cameraMex('aq_live_restart',centeredROI,binning,exposureTime)\\n')");
// 					matlabPtr->eval(u"error('centeredROI must be a 2 or 4-element vector\\n')");
// 					matlabPtr->eval(u"error('binning must be a numeric scalar with value 1, 2, 4\\n')");
// 					matlabPtr->eval(u"error('exposureTime must be a numeric scalar\\n')");
// 					return;
// 				}
                
                isROICentered = !(2 - inputs[1].getNumberOfElements());
                
                TypedArray<double> input1 = std::move(inputs[1]);
				inputROI = getPointer(input1);
                
                TypedArray<double> input2 = std::move(inputs[2]);
				binning = input2[0];
                
                TypedArray<double> input3 = std::move(inputs[3]);
				expTime = input3[0]*1e+3;
                
				switch ((int)binning) {
					case 1:
						aq_sync_prepare(isROICentered, inputROI, 1, expTime);
						break;
					case 2:
						aq_sync_prepare(isROICentered, inputROI, 2, expTime);
						break;
					case 4:
						aq_sync_prepare(isROICentered, inputROI, 4, expTime);
						break;
					default:
						matlabPtr->eval(u"error('binning must be a numeric scalar with value 1, 2, or 4\\n')");
						return;
						break;
				}
				outputs[0] = factory.createScalar<double>(exposureTimeSeconds);
				outputs[1] = factory.createScalar<double>(readoutTimeSeconds);
		// 		if( (nrhs!=3) || 
		// 			(!mxIsNumeric(inputs[1]) || mxGetM(inputs[1])!=1 || mxGetN(inputs[1])!=2 || mxGetNumberOfDimensions(inputs[1])!=2) || // 1st parameter, ROI
		// 			(!mxIsNumeric(inputs[2]) || mxGetM(inputs[2])!=1 || mxGetN(inputs[2])!=1 || mxGetNumberOfDimensions(inputs[2])!=2) ) { // 2nd parameter, exposure time
		// 			matlabPtr->eval(u"error('aq_sync_prepare requires two parameters: cameraMex('aq_sync_prepare',centeredROI,exposureTime)\n");
		// 			matlabPtr->eval(u"error('centeredROI must be a 1x2 vector and exposureTime must be a numeric scalar\n");
		// 			return;
		// 		}
		// 		inputROI = ((double *)(mxGetData(inputs[1])));
		// 		expTime = *((double *)(mxGetData(inputs[2])));
		// 		aq_sync_prepare(inputROI,expTime);
			} else if (!std::strcmp(input_command.c_str(), "aq_sync_start")) { 		////  aq_sync_start
				matlabPtr->eval(u"fprintf('running aq_sync_start\\n')");
// 				if( (inputs.size() != 8) ||
// 					(inputs[1].getDimensions()[0]!=1 || 
// 					inputs[1].getDimensions()[1] != 1 || getNumElements(inputs[1].getDimensions())!=2) ) {
// 					matlabPtr->eval(u"error('aq_sync_start requires three parameters: cameraMex('aq_sync_start',recordFrames,\\n')");
// 					matlabPtr->eval(u"error('    dcimg_fpath,dcimg_fext,\\n')");
// 					matlabPtr->eval(u"error('    bin_temp_fpath,bin_temp_fext,\\n')");
// 					matlabPtr->eval(u"error('    bin_tgt_fpath,bin_tgt_fext)\\n')");
// 					matlabPtr->eval(u"error('recordFrames must be a numeric scalar\\n')");
// 					matlabPtr->eval(u"error('fpath must be double back-slashed path strings\\n')");
// 					matlabPtr->eval(u"error('fext must be extension strings\\n')");
// 					return;
// 				}
                
                TypedArray<double> input1 = std::move(inputs[1]);
				double recFrames = input1[0];
				
                param1 = CharArray(inputs[2]).toAscii().c_str();
                param2 = CharArray(inputs[3]).toAscii().c_str();
                param3 = CharArray(inputs[4]).toAscii().c_str();
                param4 = CharArray(inputs[5]).toAscii().c_str();
                param5 = CharArray(inputs[6]).toAscii().c_str();
                param6 = CharArray(inputs[7]).toAscii().c_str();
                
				strcpy(fpath_bin_temp_copy,param3);
				strcpy(fext_bin_temp_copy,param4);
				strcpy(fpath_bin_tgt_copy,param5);
				strcpy(fext_bin_tgt_copy,param6);
				aq_sync_start(recFrames,fpath_bin_tgt_copy,"bin");
//                 aq_sync_start(recFrames,param5,param6);
				matlabPtr->eval(u"fprintf('aq_sync_start ok\\n')");
                matlabPtr->feval(u"fprintf",inputs[6]);
                matlabPtr->feval(u"fprintf",inputs[7]);
                matlabPtr->eval(u"fprintf('\\n')");
			} else if (!std::strcmp(input_command.c_str(), "test")) { 				////  test
				matlabPtr->eval(u"fprintf('testing nothing\\n')");
			} 
//             else if (!strcmp(input_command, "aq_sync_stop")) { 		////  aq_sync_stop
// 				matlabPtr->eval(u"fprintf('running aq_sync_stop\n')");
// 				aq_sync_stop();
// 			} 
            else if (!std::strcmp(input_command.c_str(), "aq_live_restart")) { 	////  aq_live_restart
				matlabPtr->eval(u"fprintf('running aq_live_restart\\n')");
// 				if( (inputs.size() != 4) || 
// 					((!inputs[1].getNumberOfElements() == 2 || inputs[1].getNumberOfElements() == 4) || getNumElements(inputs[1].getDimensions()) != 2) || // 1st parameter, ROI
// 					(inputs[2].getNumberOfElements() != 1 || getNumElements(inputs[2].getDimensions()) != 2) || // 2nd parameter, binning
// 					(inputs[3].getNumberOfElements() !=1 || getNumElements(inputs[3].getDimensions()) != 2) ) { // 3nd parameter, exposure time
// 					matlabPtr->eval(u"error('aq_live_restart requires two parameters: cameraMex(''aq_live_restart'',centeredROI,binning,exposureTime)\\n')");
// 					matlabPtr->eval(u"error('centeredROI must be a 2 or 4-element vector\\n')");
// 					matlabPtr->eval(u"error('binning must be a numeric scalar with value 1, 2, 4\\n')");
// 					matlabPtr->eval(u"error('exposureTime must be a numeric scalar\\n')");
// 					return;
// 				}
                isROICentered = !(2 - inputs[1].getNumberOfElements());
                
                TypedArray<double> input1 = std::move(inputs[1]);
				inputROI = getPointer(input1);
                
                TypedArray<double> input2 = std::move(inputs[2]);
				binning = input2[0];
                
                TypedArray<double> input3 = std::move(inputs[3]);
				expTime = input3[0]*1e+3;

				switch ((int)binning) {
					case 1:
                        matlabPtr->eval(u"fprintf('starting aq live restart \\n')");
						if(aq_live_restart(isROICentered, inputROI, 1, expTime)){
                            matlabPtr->eval(u"fprintf('aq live restart...ok \\n')");
                        }
                        else{
                            matlabPtr->eval(u"fprintf('aq live restart...failed \\n')");
                        }
                        
						break;
					case 2:
						aq_live_restart(isROICentered, inputROI, 2, expTime);
						break;
					case 4:
						aq_live_restart(isROICentered, inputROI, 4, expTime);
						break;
					default:
						matlabPtr->eval(u"error('binning must be a numeric scalar with value 1, 2, or 4\\n')");
						return;
						break;
				}

				outputs[0] = factory.createScalar<double>(exposureTimeSeconds);
				outputs[1] = factory.createScalar<double>(readoutTimeSeconds);
			} else if (!std::strcmp(input_command.c_str(), "aq_snap")) { 			////  aq_snap
				matlabPtr->eval(u"fprintf('running aq_snap\\n')");
				if (inputs.size() > 1){
					matlabPtr->eval(u"error('aq_snap requires at most one output: img = cameraMex('aq_snap')\\n')");
				}
				newFrame = aq_snap();
				if (newFrame) {
// 					WaitForSingleObject(ghMutexFrame, INFINITE);
                    ghMutexFrame.lock();
						memcpy(lastImageData, newFrame->buf, (newFrame->width) * (newFrame->height)*sizeof(uint16_t));
						lastImageDataWidth  = newFrame->width;
						lastImageDataHeight = newFrame->height;
// 					ReleaseMutex(ghMutexFrame);
                    ghMutexFrame.unlock();
					matlabPtr->eval(u"fprintf('got new image, returning\\n')");
				}else{
					matlabPtr->eval(u"fprintf('timeout with no new image, returning last\\n')");
				}

				// there is a new image
//                 TypedArray<uint16_t> outBuf = factory.createArray<uint16_t>({lastImageDataHeight, lastImageDataWidth});
                TypedArray<uint16_t> outBuf = factory.createArray<uint16_t>({lastImageDataWidth, lastImageDataHeight});
                memcpy(getPointer(outBuf), lastImageData, lastImageDataWidth*lastImageDataHeight*sizeof(uint16_t));
				
				// sprintf(msg,"dimensions: %d %d\n", lastImageDataWidth, lastImageDataHeight);
				// matlabPtr->eval(u"fprintf('msg);
				outputs[0] = outBuf;

			} else if (!std::strcmp(input_command.c_str(), "aq_live_stop")) { 		////  aq_live_stop
				matlabPtr->eval(u"fprintf('running aq_live_stop\\n')");
				aq_live_stop();
			} else if (!std::strcmp(input_command.c_str(), "verbose")) { 		////  aq_live_stop
				verboseFlag = !verboseFlag;
                matlabPtr->eval(u"fprintf('verbose flag set')");
			} else { 													////  unknown
				matlabPtr->feval(u"fprintf",inputs[0]);
			}
			return;
    }
};