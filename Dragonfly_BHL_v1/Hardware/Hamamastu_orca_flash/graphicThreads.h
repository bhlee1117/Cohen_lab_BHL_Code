/* * * * * * * * * * * * * * * * * * *\
 *                                    *
 *   2018 Vicente Parot               *
 *   Cohen Lab - Harvard University   *
 *                                    *
 \* * * * * * * * * * * * * * * * * * */

#ifndef _include_graphicThreads_h_
#define _include_graphicThreads_h_

#include <SDL.h> // graphics layer
#include "mex.h" // matlab interface
#include <windows.h> // for threads, mutex
#include <process.h> // for threads, mutex
#include "dcamapi4.h"

#include "camControl.h"

// rendering interpolation
// "0" nearest neigbor
// "1" linear filtering (OpenGL and Direct3D)
// "2" anisotropic filtering (Direct3D)
#define RENDER_SCALE_QUALITY ("0")

// rendering dimensions
#define DEFAULT_MONITOR_INDEX 1
#define DEFAULT_SCREEN_WIDTH 955
#define DEFAULT_SCREEN_HEIGHT 1174
#define DEFAULT_SCREEN_X_POSITION 2
#define DEFAULT_SCREEN_Y_POSITION 24
#define NPIX_SAT_WARNING 4

// buffers dimensions
#define IMG_WIDTH 2048
#define IMG_HEIGHT 2048
#define IMG_BYTES_PER_PIXEL 2
#define HISTO_WIDTH 955
#define HISTO_HEIGHT 256

// bit masks
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
#define rmask 0xff000000
#define gmask 0x00ff0000
#define bmask 0x0000ff00
#define amask 0x000000ff
#else                   
#define rmask 0x000000ff
#define gmask 0x0000ff00
#define bmask 0x00ff0000
#define amask 0xff000000
#endif

// // rendering areas
// static SDL_Rect imgSrcRect;
// static SDL_Rect imgDestRect;
// static SDL_Rect histoSrcRect;
// static SDL_Rect histoDestRect;

// // pointers to always-valid surfaces
// static SDL_Surface* myReadyImgSurf = NULL;
// static SDL_Surface* myReadyHistoSurf = NULL;

// // pointers to allocated surface buffers
// static SDL_Surface* myWork1ImgSurf = NULL;
// static SDL_Surface* myWork2ImgSurf = NULL;
// static SDL_Surface* myWork1HistoSurf = NULL;
// static SDL_Surface* myWork2HistoSurf = NULL;

// // rendering pointers
// static SDL_Window* window;
// static SDL_Renderer* renderer;
// static SDL_Texture* myImgTex = NULL;
// static SDL_Texture* myHistoTex = NULL;

// // pointers to pixel data buffers
// static uint32_t *workImgSurfBuffer = NULL;
// static uint32_t *workHistoSurfBuffer = NULL;
// static uint32_t *tempImgBufferRGBA = NULL;
// static uint32_t *tempHistBufferRGBA = NULL;

// // synchronization flags
// static bool continueRunningAllThreads;
// static bool imgToggleFlag = true;
// static bool histoToggleFlag = true;
// static bool isThereANewImage = true;
// static bool isThereANewCoordinate = true;

// // thread handles
// static HANDLE hImgThread;
// static HANDLE hProcThread;
// static HANDLE hHistoThread;

// // mutex handles
// static HANDLE ghMutexImg; 
// static HANDLE ghMutexHisto; 

// // pointer to test image data
// static uint16_t *tempImgBufferU16 = NULL;

// // mouse coordinates
// static int latest_mouse_x = 0;
// static int latest_mouse_y = 0;

void calc_and_update_histo(void);
void calc_and_update_image(void);
void update_rendering(void);
unsigned __stdcall imgThreadFcn  ( void* pArguments );
unsigned __stdcall procThreadFcn ( void* pArguments );
unsigned __stdcall histoThreadFcn( void* pArguments );
unsigned __stdcall convertThreadFcn( void* pArguments );
void start_SDL_window(uint16_t *windowLocationParams);
void load_test_image(void);
void launch_threads(void);
void close_SDL_window_and_threads(void);

#endif // define _include_graphicThreads_h_
