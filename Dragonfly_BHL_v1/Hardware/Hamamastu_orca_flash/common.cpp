/* * * * * * * * * * * * * * * * * * * *\
 *                                     *
 *   2018 Vicente Parot                *
 *   Cohen Lab - Harvard University    *
 *                                     *
 *   Modified from Hamamatsu DCAMAPI   *
 *                                     *
\* * * * * * * * * * * * * * * * * * * */

// console/misc/common.cpp
//

#include	"console.h"
#include	"common.h"
#include 	"mex.h"

#ifndef ASSERT
#define	ASSERT(c)
#endif

extern bool verboseFlag;

void message( const char* fmt){
	if (verboseFlag){		
		if( fmt != NULL ){
			printf( fmt );
		}
	}
}

void error( const char* fmt){
	if( fmt != NULL ){
		printf("!!! ");
		printf( fmt );
	}
}

// ----------------------------------------------------------------

inline const int my_dcamdev_string( DCAMERR& err, HDCAM hdcam, int32 idStr, char* text, int32 textbytes )
{
	DCAMDEV_STRING	param;
	memset( &param, 0, sizeof(param) );
	param.size		= sizeof(param);
	param.text		= text;
	param.textbytes	= textbytes;
	param.iString	= idStr;
	
	err = dcamdev_getstring( hdcam, &param );
	return ! failed( err );
}

// ----------------------------------------------------------------
void dcamcon_show_dcamerr( HDCAM hdcam, DCAMERR errid, const char* apiname)
{
	char msg[MAX_MSG_LEN];
	char errtext[ 256 ];

	DCAMERR err;
	my_dcamdev_string( err, hdcam, errid, errtext, sizeof(errtext) );

	sprintf(msg,"FAILED: (DCAMERR)0x%08X %s @ %s", errid, errtext, apiname );
	message(msg);

	// if( fmt != NULL )
	// {
		// printf( " : " );

		// va_list	arg;
		// va_start(arg,fmt);
		// printf( fmt, arg );
		// va_end(arg);
	// }

	printf( "\n" );
}

// void dcamcon_show_dcamerr( HDCAM hdcam, DCAMERR errid, const char* apiname)
// {
	// char errtext[ 256 ];

	// DCAMERR err;
	// my_dcamdev_string( err, hdcam, errid, errtext, sizeof(errtext) );

	// printf( "FAILED: (DCAMERR)0x%08X %s @ %s", errid, errtext, apiname );

	// printf( "\n" );
// }

// ----------------------------------------------------------------

void dcamcon_show_dcamdev_info( HDCAM hdcam )
{
	char	model[ 256 ];
	char	cameraid[ 64 ];
	char	bus[ 64 ];

	DCAMERR	err;
	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_MODEL,    model,    sizeof(model)   ) )
	{
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_MODEL)\n" );
	}
	else
	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_CAMERAID, cameraid, sizeof(cameraid)) )
	{
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_CAMERAID)\n" );
	}
	else
	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_BUS,      bus,      sizeof(bus)     ) )
	{
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_BUS)\n" );
	}
	else
	{
		printf( "%s (%s) on %s\n", model, cameraid, bus );
	}
}

// show HDCAM camera information by text.
void dcamcon_show_dcamdev_info_detail( HDCAM hdcam )
{
	char	buf[ 256 ];

	DCAMERR	err;
	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_VENDOR, buf, sizeof(buf) ) )
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_VENDOR)\n" );
	else
		printf( "DCAM_IDSTR_VENDOR         = %s\n", buf );

	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_MODEL, buf, sizeof(buf) ) )
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_MODEL)\n" );
	else
		printf( "DCAM_IDSTR_MODEL          = %s\n", buf );

	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_CAMERAID, buf, sizeof(buf) ) )
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_CAMERAID)\n" );
	else
		printf( "DCAM_IDSTR_CAMERAID       = %s\n", buf );

	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_BUS, buf, sizeof(buf) ) )
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_BUS)\n" );
	else
		printf( "DCAM_IDSTR_BUS            = %s\n", buf );


	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_CAMERAVERSION, buf, sizeof(buf) ) )
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_CAMERAVERSION)\n" );
	else
		printf( "DCAM_IDSTR_CAMERAVERSION  = %s\n", buf );

	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_DRIVERVERSION, buf, sizeof(buf) ) )
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_DRIVERVERSION)\n" );
	else
		printf( "DCAM_IDSTR_DRIVERVERSION  = %s\n", buf );

	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_MODULEVERSION, buf, sizeof(buf) ) )
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_MODULEVERSION)\n" );
	else
		printf( "DCAM_IDSTR_MODULEVERSION  = %s\n", buf );

	if( ! my_dcamdev_string( err, hdcam, DCAM_IDSTR_DCAMAPIVERSION, buf, sizeof(buf) ) )
		dcamcon_show_dcamerr( hdcam, err, "dcamdev_getstring(DCAM_IDSTR_DCAMAPIVERSION)\n" );
	else
		printf( "DCAM_IDSTR_DCAMAPIVERSION = %s\n", buf );
}


// ----------------------------------------------------------------
// initialize DCAM-API and get HDCAM camera handle.

HDCAM dcamcon_init_open()
{
	char msg[MAX_MSG_LEN];
	// message("dcamcon_init_open entry\n");
	// Initialize DCAM-API ver 4.0
	DCAMAPI_INIT	apiinit;
	memset( &apiinit, 0, sizeof(apiinit) );
	apiinit.size	= sizeof(apiinit);
	
	DCAMERR	err;
	err = dcamapi_init( &apiinit );
	if( failed( err ) ){
		// failure
		dcamcon_show_dcamerr( NULL, err, "dcamapi_init()" );
		return NULL;
	}else{
		sprintf(msg,"dcamapi_init ok\n");
		message(msg);
		int32	nDevice = apiinit.iDeviceCount;
		// message("found %d device%s\n",nDevice,(nDevice-1)?"s":"");
		if (nDevice>0){
			ASSERT( nDevice > 0 );	// nDevice must be larger than 0
			
			int32	iDevice;
			
			// show all camera information by text
			for( iDevice = 0; iDevice < nDevice; iDevice++ )
			{
				dcamcon_show_dcamdev_info( (HDCAM)(intptr_t)iDevice );
			}
			
			if( nDevice > 1 )
			{
				sprintf(msg,"selecting first camera available from the list\n",apiinit.iDeviceCount);
				message(msg);
				// choose one camera from the list if there are two or more cameras.
				// printf( "choose one of camera from above list by index (0-%d) >", nDevice-1 );
				iDevice = 0;
				// iDevice = -1;

				// char	buf[256];
				// while( fgets( buf, sizeof(buf), stdin ) != NULL )
				// {
					// if( _stricmp( buf, "exit" ) == 0 || _stricmp( buf, "exit\n" ) == 0 )
					// {
						// iDevice = -1;
						// break;
					// }

					// iDevice = atoi( buf );
					// if( 0 <= iDevice && iDevice < nDevice )
						// break;
				// }
			}
			else
			{
				iDevice = 0;
			}
			
			if( 0 <= iDevice && iDevice < nDevice )
			{
				// open specified camera
				DCAMDEV_OPEN	devopen;
				memset( &devopen, 0, sizeof(devopen) );
				devopen.size	= sizeof(devopen);
				devopen.index	= iDevice;
				err = dcamdev_open( &devopen );
				if( ! failed(err) )
				{
					message("dcamdev_open ok\n");
					
					HDCAM	hdcam = devopen.hdcam;
					
					dcamcon_show_dcamdev_info_detail( hdcam );

					// success
					return hdcam;
				}
				
				dcamcon_show_dcamerr( (HDCAM)(intptr_t)iDevice, err, "dcamdev_open()");
				printf("index is %d\n", iDevice );
			}
			
			// uninitialize DCAM-API
			dcamapi_uninit();
		}else{ // (nDevice>0)
			message("no camera found\n");
			return NULL;
		}
	}		
	

	
	// failure
	return NULL;
}
