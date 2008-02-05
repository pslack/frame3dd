/** @FILE
	parser for microstran properties.txt files
*/

#ifndef SECTIONSPARSER_H
#define SECTIONSPARSER_H

#include "config.h"
#include "sections.h"
#include "parse.h"

#ifdef __cplusplus
extern "C"{
#endif

MSTRANP_API cbool parseSections(parse *p, section_library *a);

#ifdef __cplusplus
}
#endif

#endif

