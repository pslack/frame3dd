/** @FILE
	parser for microstran .arc files
*/

#ifndef MODELPARSER_H
#define MODELPARSER_H

#include "config.h"
#include "model.h"
#include "parse.h"

#ifdef __cplusplus
extern "C"{
#endif

MSTRANP_API cbool parseModelMicrostran(parse *p, model **m);

#ifdef __cplusplus
}
#endif

#endif

