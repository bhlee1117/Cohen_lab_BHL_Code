
#ifndef pvcam_graphicThreads_h_
#define pvcam_graphicThreads_h_

#include <SDL.h> // graphics layer
// #include <windows.h> // for threads, mutex
// #include <process.h> // for threads, mutex
#include <mutex> // for threads, mutex
#include <thread>
#include <chrono>

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
#define IMG_WIDTH 3200
#define IMG_HEIGHT 3200
#define IMG_BYTES_PER_PIXEL 2
#define HISTO_WIDTH 955
#define HISTO_HEIGHT 256
#define HISTO_SIZE 512

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


void calc_and_update_histo(void);
void calc_and_update_image(void);
void update_rendering(void);
void imgThreadFcn();
void procThreadFcn();
void histoThreadFcn();
void start_SDL_window(uint16_t *windowLocationParams);
void load_test_image(void);
void launch_threads(void);
bool close_SDL_window_and_threads(void);


#endif // define _include_graphicThreads_h_
