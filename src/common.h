#ifndef FRAME_COMMON_H
#define FRAME_COMMON_H

/* this file contains some #defines to set up 'float' to be 'double' instead. */

#ifndef PI
#define PI	3.141592653589793
#endif

#define Zvert 1		/* Zvert=0: Y axis is vertical ... rot Z, then rot Y */
			/* Zvert=1: Z axis is vertical ... rot Y, then rot Z */

/* ----------------------- double precision --------------------------------- */
#define float double
#define vector dvector
#define matrix dmatrix
#define free_vector free_dvector
#define free_matrix free_dmatrix
#define save_matrix save_dmatrix
#define show_vector show_dvector 
#define save_ut_matrix save_ut_dmatrix
/* ----------------------- double precision --------------------------------- */
/* also modify all fscanf lines ...     :1,$ s/%f/%lf/g    :1,$ s/%lf/%f/g    */

#endif /* FRAME_COMMON_H */

