/* * * * * * * * * * * * * * * * * * *\
	Header file for cameraMex.cpp
	cameraMex.cpp: matlab wrapper functions
	for PVCAM and live image display

	2022 Yitong Qi
	Cohen Lab - Harvard University
 \* * * * * * * * * * * * * * * * * * */

#ifndef cameraMex_hpp_
#define cameraMex_hpp_

// system library
#include <cstdlib>
#include <utility>
#include <cstring>
// #include <windows.h>
// #include <process.h>
#include <mutex>
#include <thread>
#include <chrono>
#include <cstdio>

// matlab library
#include "mex.hpp"
#include "mexAdapter.hpp"

// custom library
#include "Common.h"
#include "pvcamControl.h"
#include "pvcam_graphicThreads.h"


#endif // define _include_cameraMex_h_