/* error.c */

#include "types.h"
#include "print.h"
#include "error.h"

/* global variables */

print *errorGlobalPrint = NULL;

/* functions */

void errorInitPrint( void )
{
  if ( !errorGlobalPrint )
    errorGlobalPrint = printCreateFile( stderr );
}

void errorSetPrint( print *p, cbool dispose )
{
  if ( dispose && errorGlobalPrint )
    printDispose( errorGlobalPrint );
  errorGlobalPrint = p;
}

print *errorGetPrint( void )
{
  return errorGlobalPrint;
}

/* end error.c */
