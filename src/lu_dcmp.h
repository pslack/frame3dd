#ifndef FRAME_LU_DCMP_H
#define FRAME_LU_DCMP_H

/**
	Solves [A]{x} = {b}, simply and efficiently, by performing an 
	LU-decomposition of matrix [A]. No pivoting is performed.
	@param A is a diagonally dominant matrix of dimension [1..n][1..n]. 
	@param b is a r.h.s. vector of dimension [1..n].

	{b} is updated using [LU] and then back-substitution is done to obtain {x}.  
	{b} is replaced by {x} and [A] is replaced by the LU-reduction of itself.
*/
void lu_dcmp (
		float **A /* the system matrix, and its LU-reduction */
		, int n /* the dimension of the matrix */	
		, float *b /* the right hand side vector, and the solution vector	*/
		, int reduce /* 1: do a forward reduction; 0: don't do the reduction */
		, int solve /* 1: do a back substitution for {x};  0: do no bk-sub'n */
		, int *pd /* 1: positive diagonal  and  successful LU decomp'n	*/
);

#endif /* FRAME_LU_DCMP_H */
