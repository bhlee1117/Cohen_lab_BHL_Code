/******************************************************************************/
/* Copyright (C) Teledyne Photometrics. All rights reserved.                  */
/******************************************************************************/
#pragma once
#ifndef VERSION_H_
#define VERSION_H_

////////////////////////////////////////////////////////////////////////////////
//////// Update minor/major versions with more significant changes /////////////
#define VERSION_MAJOR 1
#define VERSION_MINOR 6
////////////////////////////////////////////////////////////////////////////////
//////// Increase build version for every change that goes to SVN //////////////
#define system_version_build 63
////////////////////////////////////////////////////////////////////////////////

// Nicer name alias (system_version_build name was required by build_count.exe)
#define VERSION_BUILD system_version_build

// Stringifying macros
#define STR_EXPAND(input) #input
#define STR(input) STR_EXPAND(input)

// Handy expand macros
#define VERSION_NUMBER VERSION_MAJOR,VERSION_MINOR,VERSION_BUILD
#define VERSION_NUMBER_STR STR(VERSION_MAJOR) "." STR(VERSION_MINOR) "." STR(VERSION_BUILD)

#define COMPANY_NAME_STR "Teledyne Photometrics"
#define LEGAL_COPYRIGHT_STR "Copyright (C) Teledyne Digital Imaging US, Inc."

#endif // VERSION_H_
