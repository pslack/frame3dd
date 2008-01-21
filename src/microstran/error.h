#ifndef ERROR_H
#define ERROR_H

#include "types.h"
#include "print.h"

/* global variables */

extern print *errorGlobalPrint;

/* macros */

#define errorReport(args)        errorReportLocationWith(1,__FILE__, __LINE__,args)
#define errorReportWith(ex,args) errorReportLocationWith(ex,__FILE__, __LINE__,args)
#define errorReportLocationWith(ex,file,line,args) \
  ( errorInitPrint() \
  , printFormat( errorGlobalPrint, ("error in %s, line %d: ", file, line) ) \
  , printFormat( errorGlobalPrint, args ) \
  , exit( ex ) \
  )

#define errorReportExt(args)        errorReportExtWith(1,args)
#define errorReportExtWith(ex,args) \
  ( errorInitPrint() \
  , printFormat( errorGlobalPrint, ("error: ") ) \
  , printFormat( errorGlobalPrint, args ) \
  , exit( ex ) \
  )

/* functions */

void  errorInitPrint( void );
void  errorSetPrint( print *p, cbool dispose );
print *errorGetPrint( void );

#endif /* ERROR_H */

