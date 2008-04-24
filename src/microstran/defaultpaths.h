#ifndef FRAME3DD_DEFAULTPATH_H
#define FRAME3DD_DEFAULTPATH_H

#include "config.h"


#ifdef __cplusplus
extern "C"{
#endif

/**
	Get a file-system-specific location for data files associated with
	the microstran parser (specifically, 'properties.txt').

	You do not own the returned string, and you can not depend
	on it remaining unmodified. You should copy the returned
	string to somewhere safe.
*/
MSTRANP_API const char *get_default_data_path();


#ifdef __cplusplus
};
#endif

#endif

