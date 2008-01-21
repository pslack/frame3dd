/** @FILE
	parser for microstran .arc files
*/
#ifndef DISPLACEMENTPARSER_H
#define DISPLACEMENTPARSER_H

#include "config.h"
#include "displacements.h"
#include "parse.h"

#ifdef __cplusplus
extern "C"{
#endif

MSTRANP_API cbool parseMicrostranDisplacements(parse *p, modeldisplacements **m);

#ifdef __cplusplus
}
#endif

#endif
