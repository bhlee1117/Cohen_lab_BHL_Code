// Include utilities for converting mxGPUArray to cv::cuda::GpuMat
#include "matrix.h"
#include "opencvgpumex.hpp"
#include "opencvmex.hpp"
#include "opencv2/cudaarithm.hpp"
// #include "opencv2/cudafeatures2d.hpp"
#include "opencv2/cudaoptflow.hpp"
#include "gpu/mxGPUArray.h"


using namespace cv;

///////////////////////////////////////////////////////////////////////////
// Main entry point to a MEX function
///////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{  
    
    // Ensure MATLAB's GPU support is available.
    mxInitGPU();

//     checkInputs(nrhs, prhs);
    Ptr<cuda::FarnebackOpticalFlow> OF = cuda::FarnebackOpticalFlow::create();
    OF->setNumIters(3); OF->setNumLevels(3); // farneback
//     OF->setWinSize(20); OF->setFlags(OPTFLOW_FARNEBACK_GAUSSIAN); // farneback
//     Ptr<cuda::BroxOpticalFlow> OF = cuda::BroxOpticalFlow::create(20, 5, 0.8,5,3,10);
//     Ptr<cuda::OpticalFlowDual_TVL1> OF = cuda::OpticalFlowDual_TVL1::create();
//     Ptr<cuda::DensePyrLKOpticalFlow> OF = cuda::DensePyrLKOpticalFlow::create();
    
    Ptr<cv::cuda::GpuMat> gpuMov = ocvMxGpuArrayToGpuMat_single(prhs[1]); // Convert gpuArray into cv::gpu::GpuMat 
    Ptr<cv::cuda::GpuMat> gpuMov_template = ocvMxGpuArrayToGpuMat_single(prhs[0]); // Convert gpuArray into cv::gpu::GpuMat 
    cuda::GpuMat gFlow(gpuMov->size(),CV_32FC2);
    
    OF->calc(*gpuMov,*gpuMov_template,gFlow);
    
    plhs[0] = ocvMxGpuArrayFromGpuMat_single(gFlow.reshape(1));
//     gFlow.release();
}


