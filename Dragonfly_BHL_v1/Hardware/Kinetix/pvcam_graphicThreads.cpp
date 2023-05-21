/* * * * * * * * * * * * * * * * * * *\
 *                                    *
 *   2018 Vicente Parot               *
 *	 2022 Yitong Qi					  *
 *   Cohen Lab - Harvard University   *
 *                                    *
 \* * * * * * * * * * * * * * * * * * */

#include "pvcam_graphicThreads.h"

#include <cstring>
#include <stdlib.h>
#include <stdio.h>

#include "Common.h"
#include "pvcamControl.h"

// camera image access
extern bool verboseFlag;
extern bool isCapturing;
extern bool isRecording;
// extern DCAMBUF_FRAME *lastBufferFrame;
// extern DCAMBUF_FRAME *lastValidBufferFrame;
extern uint16_t *lastImageData;
extern uint16_t lastImageDataWidth;
extern uint16_t lastImageDataHeight;

// rendering areas
static SDL_Rect imgSrcRect;
static SDL_Rect imgDestRect;
static SDL_Rect histoSrcRect;
static SDL_Rect histoDestRect;

// pointers to always-valid surfaces
static SDL_Surface* myReadyImgSurf = NULL;
static SDL_Surface* myReadyHistoSurf = NULL;

// pointers to allocated surface buffers
static SDL_Surface* myWork1ImgSurf = NULL;
static SDL_Surface* myWork2ImgSurf = NULL;
static SDL_Surface* myWork1HistoSurf = NULL;
static SDL_Surface* myWork2HistoSurf = NULL;

// rendering pointers
static SDL_Window* window;
static SDL_Renderer* renderer;
static SDL_Texture* myImgTex = NULL;
static SDL_Texture* myHistoTex = NULL;

// pointers to pixel data buffers
static uint32_t *workImgSurfBuffer = NULL;
static uint32_t *workHistoSurfBuffer = NULL;
static uint32_t *tempImgBufferRGBA = NULL;
static uint32_t *tempHistBufferRGBA = NULL;

// synchronization flags
static bool continueRunningAllThreads;
static bool imgToggleFlag = true;
static bool histoToggleFlag = true;
static bool needsHistUpdate = true;
static bool needsImgUpdate = true;

// thread handles
static HANDLE hImgThread;
static HANDLE hProcThread;
static HANDLE hHistoThread;
HANDLE hConvertThread;

// mutex handles
HANDLE ghMutexImg; 
HANDLE ghMutexLims; 
HANDLE ghMutexHisto; 
HANDLE ghMutexFrame; 
HANDLE ghMutexCapturing; 

// image display thresholds
float updateContribution = .5; // weight of the last frame relative to the previous values
float loPctileSat = .0005; 		// fixed parameter
float hiPctileSat = .9995;
float loThresholdCounts = 0; 	// updates dynamically from histogram
float hiThresholdCounts = 65535;

uint16_t MONITOR_INDEX;
uint16_t SCREEN_WIDTH;
uint16_t SCREEN_HEIGHT;
uint16_t SCREEN_X_POSITION;
uint16_t SCREEN_Y_POSITION;

int compare( const void* a, const void* b)
{
     uint16_t int_a = * ( (uint16_t *) a );
     uint16_t int_b = * ( (uint16_t *) b );

     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

void calc_and_update_histo(void){
	int i, ix, iy, ox, oy, val;
	const int HISTO_LONG_WIDTH = USHRT_MAX+1;
	float shortHisto[HISTO_WIDTH];
	float longHisto [HISTO_LONG_WIDTH];
	
	static uint16_t mySortedSomeData[IMG_WIDTH*IMG_HEIGHT];

	static uint16_t myImageData[IMG_WIDTH*IMG_HEIGHT];
	static int myImageDataWidth = IMG_WIDTH;
	static int myImageDataHeight = IMG_HEIGHT;
	
	// read only when outside data not being written
	WaitForSingleObject(ghMutexFrame, INFINITE);
		// read last image data to local copy
		memcpy(myImageData, lastImageData, lastImageDataWidth*lastImageDataHeight*sizeof(uint16_t));
		myImageDataWidth  = lastImageDataWidth ;
		myImageDataHeight = lastImageDataHeight;
	ReleaseMutex(ghMutexFrame);
	
	// downsample xN if image size larger than (N-1)*512 to calculate saturation quickly
	int step_w = (int)(((float)myImageDataWidth -1)/(IMG_WIDTH /4) + 1);
	int step_h = (int)(((float)myImageDataHeight-1)/(IMG_HEIGHT/4) + 1);
	int myNumberOfSortedData = 0;
	for (iy = 1; iy < myImageDataHeight-1; iy++){
		if((iy%step_h)>0){
			continue;
		}
		for (ix = 0; ix < myImageDataWidth; ix++){
			if((ix%step_w)>0){
				continue;
			}
			mySortedSomeData[myNumberOfSortedData++] = myImageData[iy * myImageDataWidth + ix];
		}
	}
	// sort image values
	qsort( mySortedSomeData, myNumberOfSortedData, sizeof(uint16_t), compare );
	float di;
	// update saturation values
	WaitForSingleObject(ghMutexLims, INFINITE);
		loThresholdCounts  = loThresholdCounts *(1-updateContribution) + updateContribution*mySortedSomeData[((int)(loPctileSat*myNumberOfSortedData))%myNumberOfSortedData];
		hiThresholdCounts  = hiThresholdCounts *(1-updateContribution) + updateContribution*mySortedSomeData[((int)(hiPctileSat*myNumberOfSortedData))%myNumberOfSortedData];
	ReleaseMutex(ghMutexLims);

	// draw a screen resolution histogram 
	for(i = 0; i < HISTO_WIDTH; i++){
		shortHisto[i] = 0;
	}
	for(i = 0; i < HISTO_LONG_WIDTH; i++){
		longHisto[i] = 0;
	}
	int minh = 0; // histogram drawing limits
	int maxh = 65535;
	int tox = IMG_WIDTH/2-myImageDataWidth/2; // target offset x
	int toy = IMG_HEIGHT/2-myImageDataHeight/2; // target offset x
	for (iy = 0; iy < myImageDataHeight; iy++){
		for (ix = 0; ix < myImageDataWidth; ix++){
			val = (int)(((float)myImageData[iy * myImageDataWidth + ix]-minh)/(maxh-minh)*HISTO_WIDTH);
			if (0<=val && val<HISTO_WIDTH){
				shortHisto[val]++;
			}
			val = myImageData[iy * myImageDataWidth + ix];
			if (0<=val && val<HISTO_LONG_WIDTH){
				longHisto[val]++;
			}
		}
	}
	// normalize histogram to its max 
	float histoMax = 0;
	for(i = 1; i < HISTO_WIDTH-1; i++){
		histoMax = (shortHisto[i]>histoMax)?shortHisto[i]:histoMax;
	}
	for (ix = 0; ix < HISTO_WIDTH; ix++){
		for (iy = 0; iy < HISTO_HEIGHT; iy++){
			// mind bit-endian-ness here
			float step = (float)HISTO_HEIGHT/histoMax;
			float diff = (shortHisto[ix]*step) - (HISTO_HEIGHT-iy-1);
			diff = (  -1>diff)?-1:diff;
			diff = (diff>0   )? 0:diff;
			diff = diff*(255-32)+255;
			val = (int)diff;
			tempHistBufferRGBA[iy * HISTO_WIDTH + ix] = (0xFF << 24) | (val << 16) | (val << 8) | val;
		}
	}
	// draw vertical bars on the saturation values
	ix = (int)(loThresholdCounts/65536*HISTO_WIDTH);
	ix = ix%HISTO_WIDTH;
	for (iy = 0; iy < HISTO_HEIGHT; iy++){
		// mind bit-endian-ness here
		val = 127;
		tempHistBufferRGBA[iy * HISTO_WIDTH + ix] = (0xFF << 24) | (val << 16) | (val << 8) | val;
		tempHistBufferRGBA[iy * HISTO_WIDTH + ix+1] = tempHistBufferRGBA[iy * HISTO_WIDTH + ix];
	}
	ix = (int)(hiThresholdCounts/65536*HISTO_WIDTH); 
	ix = ix%HISTO_WIDTH;
	for (iy = 0; iy < HISTO_HEIGHT; iy++){
		// mind bit-endian-ness here
		val = 127;
		tempHistBufferRGBA[iy * HISTO_WIDTH + ix] = (0xFF << 24) | (val << 16) | (val << 8) | val;
        tempHistBufferRGBA[iy * HISTO_WIDTH + ix-1] = tempHistBufferRGBA[iy * HISTO_WIDTH + ix];
	}
	// draw vertical red bars on camera saturation
	if (longHisto [HISTO_LONG_WIDTH-1]>0) {
		for (ix = HISTO_WIDTH-NPIX_SAT_WARNING; ix < HISTO_WIDTH; ix++){ // N px thick
			for (iy = 0; iy < HISTO_HEIGHT; iy++){
				// mind bit-endian-ness here
				val = 255;
				tempHistBufferRGBA[iy * HISTO_WIDTH + ix] = (0xFF << 24) | (0 << 16) | (0 << 8) | val;
			}
		}
	}
	if (longHisto [0]>0) {
		for (ix = 0; ix < NPIX_SAT_WARNING; ix++){ // N px thick
			for (iy = 0; iy < HISTO_HEIGHT; iy++){
				// mind bit-endian-ness here
				val = 255;
				tempHistBufferRGBA[iy * HISTO_WIDTH + ix] = (0xFF << 24) | (0 << 16) | (0 << 8) | val;
			}
		}
	}
	
	histoToggleFlag = !histoToggleFlag;
	workHistoSurfBuffer = (uint32_t *)(histoToggleFlag?myWork1HistoSurf:myWork2HistoSurf)->pixels;
	memcpy( workHistoSurfBuffer, tempHistBufferRGBA, HISTO_WIDTH*HISTO_HEIGHT*sizeof(uint32_t));
	
	WaitForSingleObject(ghMutexHisto, INFINITE);
		myReadyHistoSurf = histoToggleFlag?myWork1HistoSurf:myWork2HistoSurf;
	ReleaseMutex(ghMutexHisto);
}

void calc_and_update_image(void){
	int ix, iy, ox, oy;
	
	// local data copies to lock only during copy
	static uint16_t myImageData[IMG_WIDTH*IMG_HEIGHT];
	static int myImageDataWidth = IMG_WIDTH;
	static int myImageDataHeight = IMG_HEIGHT;
	
	// ensure lastImageData read while not writing elsewhere
	WaitForSingleObject(ghMutexFrame, INFINITE);
		// make local copy of latest image
		memcpy(myImageData, lastImageData, lastImageDataWidth*lastImageDataHeight*sizeof(uint16_t));
		myImageDataWidth  = lastImageDataWidth ;
		myImageDataHeight = lastImageDataHeight;
	ReleaseMutex(ghMutexFrame);

	// toggle and select which buffer to work on
	imgToggleFlag = !imgToggleFlag;
	workImgSurfBuffer = (uint32_t *)(imgToggleFlag?myWork1ImgSurf:myWork2ImgSurf)->pixels;

	// offsets due to image resizing
	int tox = IMG_WIDTH/2-myImageDataWidth/2; // target offset x
	int toy = IMG_HEIGHT/2-myImageDataHeight/2; // target offset y

	// ensure loThresholdCounts and hiThresholdCounts read while not writing elsewhere
	WaitForSingleObject(ghMutexLims, INFINITE);
		// scale and copy image data into 24 bpp surface buffer
		for (iy = 0; iy < myImageDataHeight; iy++){
			for (ix = 0; ix < myImageDataWidth; ix++){
				// mind bit-endian-ness here
				// uint8_t val = (uint8_t)(((float)myImageData[iy * myImageDataWidth + ix])/256.0);
				float fval = (float)myImageData[iy * myImageDataWidth + ix];
				fval = (fval<loThresholdCounts)?loThresholdCounts:fval;
				fval = (fval>hiThresholdCounts)?hiThresholdCounts:fval;
				uint8_t val = (uint8_t)((fval-loThresholdCounts)/(hiThresholdCounts-loThresholdCounts)*255.999);
				workImgSurfBuffer[(iy+toy) * IMG_WIDTH + (ix+tox)] = (0xFF << 24) | (val << 16) | (val << 8) | val;
			}
		}
	ReleaseMutex(ghMutexLims);
	
	// adjust drawing areas to reflect changes in image size
	imgSrcRect.x = tox;
	imgSrcRect.y = toy;
	imgSrcRect.w = myImageDataWidth;
	imgSrcRect.h = myImageDataHeight;
	float whr = (float)myImageDataWidth/myImageDataHeight; // aspect ratio
	float hwr = 1./whr;
	if (hwr < 1){
		imgDestRect.x = 0;
		imgDestRect.y = SCREEN_WIDTH*(1 - hwr)/2;
		imgDestRect.w = SCREEN_WIDTH;
		imgDestRect.h = SCREEN_WIDTH*hwr;
	}else{
		imgDestRect.x = SCREEN_WIDTH*(1 - whr)/2;
		imgDestRect.y = 0;
		imgDestRect.w = SCREEN_WIDTH*whr;
		imgDestRect.h = SCREEN_WIDTH;
	}

	// ensure myReadyImgSurf write while not reading elsewhere
	WaitForSingleObject(ghMutexImg, INFINITE);
		// set surface pointer to the buffer we just worked on
		myReadyImgSurf = imgToggleFlag?myWork1ImgSurf:myWork2ImgSurf;
	ReleaseMutex(ghMutexImg);
}

void update_rendering(void){
	SDL_RenderClear(renderer); // blank the rendering
	if(myReadyImgSurf){ // in case these pointers are NULL		
		// ensure myReadyImgSurf read while not writing elsewhere
		WaitForSingleObject(ghMutexImg, INFINITE);
			myImgTex = SDL_CreateTextureFromSurface( renderer, myReadyImgSurf);
		ReleaseMutex(ghMutexImg);
		SDL_RenderCopy(renderer, myImgTex, &imgSrcRect, &imgDestRect);
		SDL_DestroyTexture(myImgTex); // free texture
	}
	if(myReadyHistoSurf){
		// ensure myReadyHistoSurf read while not writing elsewhere
		WaitForSingleObject(ghMutexHisto, INFINITE);
			myHistoTex = SDL_CreateTextureFromSurface( renderer, myReadyHistoSurf);
		ReleaseMutex(ghMutexHisto);
		SDL_RenderCopy(renderer, myHistoTex, &histoSrcRect, &histoDestRect);
		SDL_DestroyTexture(myHistoTex); // free texture
	}
	SDL_RenderPresent(renderer); // put the rendering to screen
}


unsigned __stdcall imgThreadFcn( void* pArguments ) {
	SDL_Event event;
	bool done = false;
	static BUF_FRAME *newFrame = NULL;
	static uint16_t myImageData[IMG_WIDTH*IMG_HEIGHT];
	static int myImageDataWidth = IMG_WIDTH;
	static int myImageDataHeight = IMG_HEIGHT;

	while ( continueRunningAllThreads && !done) {
		SDL_Delay(2);
		// process window events
		while(SDL_PollEvent(&event)){
			switch (event.type){
				case SDL_QUIT:
					SDL_HideWindow(window);				
					done = true;
					break;
			}
		}
		// obtain live image onto local copy
		if (isCapturing){ 
			// read outside image only when not being modified
			WaitForSingleObject(ghMutexCapturing, INFINITE);
				if (isCapturing){
					newFrame = aq_snap();
					if (newFrame) {
						memcpy(myImageData, newFrame->buf, newFrame->width*newFrame->height*sizeof(uint16_t));
						myImageDataWidth  = newFrame->width;
						myImageDataHeight = newFrame->height;
					}
				}else{
					newFrame = NULL;
				}
			ReleaseMutex(ghMutexCapturing);
		}
		// broadcast live image from local copy if updated
		if (newFrame) {
			// write outside copy only if not being read
			WaitForSingleObject(ghMutexFrame, INFINITE);
				memcpy(lastImageData, myImageData, myImageDataWidth*myImageDataHeight*sizeof(uint16_t));
				lastImageDataWidth  = myImageDataWidth;
				lastImageDataHeight = myImageDataHeight;
				needsImgUpdate = true;
				needsHistUpdate = true;
			ReleaseMutex(ghMutexFrame);
		}
		update_rendering();	// self explanatory
	}
	_endthreadex( 0 );
	return 0;
}

unsigned __stdcall procThreadFcn( void* pArguments ) {
	while ( continueRunningAllThreads ) {
		SDL_Delay(2);
		// process image
		if (needsImgUpdate){
			calc_and_update_image();
			needsImgUpdate = false;
		}
	}
	_endthreadex( 0 );
	return 0;
}

unsigned __stdcall histoThreadFcn( void* pArguments ) {
	while ( continueRunningAllThreads ) {
		SDL_Delay(2);
		// update histogram
		if (needsHistUpdate){
			calc_and_update_histo();
			needsHistUpdate = false;
		}
	}
	_endthreadex( 0 );
	return 0;
}


void start_SDL_window(uint16_t *windowLocationParams){
		int ix, iy;

		/* Initialize SDL. */
		if (SDL_Init(SDL_INIT_VIDEO) < 0)
			return;

		// allocate img buffer
		tempImgBufferRGBA = (uint32_t *)malloc(IMG_WIDTH*IMG_HEIGHT*sizeof(uint32_t));
		tempHistBufferRGBA = (uint32_t *)malloc(HISTO_WIDTH*HISTO_HEIGHT*sizeof(uint32_t));
		
		if(windowLocationParams){
			MONITOR_INDEX     = windowLocationParams[0];
            SCREEN_WIDTH      = windowLocationParams[1];
            SCREEN_HEIGHT     = windowLocationParams[2];
            SCREEN_X_POSITION = windowLocationParams[3];
            SCREEN_Y_POSITION = windowLocationParams[4];
		}else{
			MONITOR_INDEX     = DEFAULT_MONITOR_INDEX    ;
			SCREEN_WIDTH      = DEFAULT_SCREEN_WIDTH     ;
			SCREEN_HEIGHT     = DEFAULT_SCREEN_HEIGHT    ;
			SCREEN_X_POSITION = DEFAULT_SCREEN_X_POSITION;
			SCREEN_Y_POSITION = DEFAULT_SCREEN_Y_POSITION;
		}
		SDL_Rect displayBounds;
		SDL_GetDisplayBounds( MONITOR_INDEX, &displayBounds );
		printf("window position set to %d %d %d %d on monitor %d\n",
			displayBounds.x + SCREEN_X_POSITION,
			displayBounds.y + SCREEN_Y_POSITION,
			SCREEN_WIDTH,
			SCREEN_HEIGHT,
			MONITOR_INDEX);
		
		/* Create the window where we will draw. */
		window = SDL_CreateWindow( "Dragonfly live image",
			displayBounds.x + SCREEN_X_POSITION, 
			displayBounds.y + SCREEN_Y_POSITION, 
			SCREEN_WIDTH, 
			SCREEN_HEIGHT, 
			SDL_WINDOW_SHOWN );

		/* We must call SDL_CreateRenderer in order for draw calls to affect this window. */
		renderer = SDL_CreateRenderer(window, -1, 0);

		SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, RENDER_SCALE_QUALITY);

		// Set Canvas background to dark gray 
		SDL_SetRenderDrawColor(renderer, 64, 64, 64, 255);

		myWork1ImgSurf = SDL_CreateRGBSurface(0, IMG_WIDTH, IMG_HEIGHT, 32, rmask, gmask, bmask, amask);
		if (myWork1ImgSurf == NULL){
			SDL_Log("SDL_CreateRGBSurface() failed: %s", SDL_GetError());
			exit(1);
		}
		myWork2ImgSurf = SDL_CreateRGBSurface(0, IMG_WIDTH, IMG_HEIGHT, 32, rmask, gmask, bmask, amask);
		if (myWork2ImgSurf == NULL){
			SDL_Log("SDL_CreateRGBSurface() failed: %s", SDL_GetError());
			exit(1);
		}
		
		myWork1HistoSurf = SDL_CreateRGBSurface(0, HISTO_WIDTH, HISTO_HEIGHT, 32, rmask, gmask, bmask, amask);
		if (myWork1HistoSurf == NULL){
			SDL_Log("SDL_CreateRGBSurface() failed: %s", SDL_GetError());
			exit(1);
		}
		myWork2HistoSurf = SDL_CreateRGBSurface(0, HISTO_WIDTH, HISTO_HEIGHT, 32, rmask, gmask, bmask, amask);
		if (myWork2HistoSurf == NULL){
			SDL_Log("SDL_CreateRGBSurface() failed: %s", SDL_GetError());
			exit(1);
		}
		
		imgSrcRect.x = 0;
		imgSrcRect.y = 0;
		imgSrcRect.w = IMG_WIDTH;
		imgSrcRect.h = IMG_HEIGHT;

		imgDestRect.x = 0;
		imgDestRect.y = 0;
		imgDestRect.w = SCREEN_WIDTH;
		imgDestRect.h = SCREEN_WIDTH; // make a square inside of the rectangular screen

		histoSrcRect.x = 0;
		histoSrcRect.y = 0;
		histoSrcRect.w = HISTO_WIDTH;
		histoSrcRect.h = HISTO_HEIGHT;

		histoDestRect.x = 0;
		histoDestRect.y = SCREEN_WIDTH+1;
		histoDestRect.w = SCREEN_WIDTH;
		histoDestRect.h = SCREEN_HEIGHT-SCREEN_WIDTH;
		
		// Copy image into RGB format for display
		// store a reference copy in temp buffer
		for (iy = 0; iy < lastImageDataHeight; iy++){
			for (ix = 0; ix < lastImageDataWidth; ix++){
				// mind bit-endian-ness here
				uint8_t val = (uint8_t)(((float)lastImageData[iy * lastImageDataWidth + ix])/256.0);
				tempImgBufferRGBA[iy * lastImageDataWidth + ix] = (0xFF << 24) | (val << 16) | (val << 8) | val;
			}
		}
		
		// here to try once and debug with printfs, but they will be run from the threads
		calc_and_update_histo();
		mexPrintf( "histo ok\n" );
		calc_and_update_image();	
		mexPrintf( "display ok\n" );
}

void load_test_image(void){
		// static DCAMBUF_FRAME	bufframe;
		// memset( &bufframe, 0, sizeof(bufframe) );
		// bufframe.size		= sizeof(bufframe);
		// bufframe.iFrame		= -1;				// latest captured image

		// allocate img buffer
		lastImageData = (uint16_t *)malloc(IMG_WIDTH*IMG_HEIGHT*sizeof(uint16_t));
		for (int iy = 0; iy < IMG_HEIGHT; iy++){
			for (int ix = 0; ix < IMG_WIDTH; ix++){
				lastImageData[iy * IMG_WIDTH + ix] = 16384;
			}
		}
		lastImageDataWidth  = IMG_WIDTH ;
		lastImageDataHeight = IMG_HEIGHT;
		// bufframe.buf = tempImgBufferU16;
		// bufframe.width = IMG_WIDTH;
		// bufframe.height = IMG_HEIGHT;
		// bufframe.rowbytes = IMG_WIDTH*2;
		// bufframe.type = DCAM_PIXELTYPE_MONO16;
		
		// return &bufframe;
		
		// FILE *fp;
		// // Read test image
		// fp = fopen("testimg.raw","rb");
		// fread(tempImgBufferU16, IMG_BYTES_PER_PIXEL, IMG_WIDTH*IMG_HEIGHT, fp);
		// printf("read %d %d\n",tempImgBufferU16[0],tempImgBufferU16[14]);
		// fclose(fp);
}

void launch_threads(){
	unsigned threadID;
	ghMutexImg = CreateMutex(NULL, FALSE, NULL);
	ghMutexLims = CreateMutex(NULL, FALSE, NULL);
	ghMutexHisto = CreateMutex(NULL, FALSE, NULL);
	ghMutexFrame = CreateMutex(NULL, FALSE, NULL);
	ghMutexCapturing = CreateMutex(NULL, FALSE, NULL);

	/* Create the second thread. */
	mexPrintf( "Creating second thread...\n" );
	continueRunningAllThreads = true; // make sure threads will run even after start-stop cycle
	hImgThread   = (HANDLE)_beginthreadex( NULL, 0, &imgThreadFcn, NULL, 0, &threadID );
	hProcThread  = (HANDLE)_beginthreadex( NULL, 0, &procThreadFcn, NULL, 0, &threadID );
	hHistoThread = (HANDLE)_beginthreadex( NULL, 0, &histoThreadFcn, NULL, 0, &threadID );
	mexPrintf( "beginthread called ok\n" );
}

void close_SDL_window_and_threads(){
	// hide window first, avoids delays and deadlock
	SDL_HideWindow(window);

	// stop thread to finish access to graphics components
	continueRunningAllThreads = false; // gently stop threads
	WaitForSingleObject( hImgThread  , INFINITE );
	WaitForSingleObject( hProcThread , INFINITE );
	WaitForSingleObject( hHistoThread, INFINITE );
	WaitForSingleObject( hConvertThread, INFINITE );
	CloseHandle( hImgThread   );
	CloseHandle( hProcThread   );
	CloseHandle( hHistoThread );
	CloseHandle( hConvertThread );
	CloseHandle(ghMutexImg);
	CloseHandle(ghMutexLims);
	CloseHandle(ghMutexHisto);
	CloseHandle(ghMutexFrame);
	CloseHandle(ghMutexCapturing);

	// then terminate graphics
	SDL_FreeSurface(myWork1ImgSurf  ); myWork1ImgSurf   = NULL;
	SDL_FreeSurface(myWork2ImgSurf  ); myWork2ImgSurf   = NULL;
	SDL_FreeSurface(myWork1HistoSurf); myWork1HistoSurf = NULL;
	SDL_FreeSurface(myWork2HistoSurf); myWork2HistoSurf = NULL;
	SDL_Quit();

	// free allocated memory
	// free(tempImgBufferU16);
	free(tempImgBufferRGBA);
	// tempImgBufferU16 = NULL;
	tempImgBufferRGBA = NULL;
}