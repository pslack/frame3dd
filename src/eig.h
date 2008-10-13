/*	FRAME3DD: Static and dynamic structural analysis of 2D & 3D frames and trusses
 Copyright (C) 1992-2008  Henri P. Gavin
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*//** @file
	Routines to solve the generalized eigenvalue problem

	H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu  1 March 2007
	Bathe, Finite Element Procecures in Engineering Analysis, Prentice Hall, 1982
*/
#ifndef FRAME_EIG_H
#define FRAME_EIG_H

/**
	Find the lowest m eigenvalues, w, and eigenvectors, V, of the 
	general eigenproblem, K V = w M V, using sub-space / Jacobi iteration.

	@param K is an n by n  symmetric real (stiffness) matrix
	@param M is an n by n  symmetric positive definate real (mass) matrix
	@param w is a diagonal matrix of eigen-values
	@param V is a  rectangular matrix of eigen-vectors
*/
void subspace(
	float **K, float **M
	, int n, int m /**< DoF and number of required modes	*/
	, float *w, float **V
	, float tol, float shift
	, int *iter /**< sub-space iterations */
	, int *ok /**< Sturm check result */
);


/**
	carry out matrix-matrix-matrix multiplication for symmetric A
	C = X' A X     C is J by J	X is N by J	A is N by N
*/
void xtAx(float **A, float **X, float **C, int N, int J);

/**
	calculate the lowest m eigen-values and eigen-vectors of the
	generalized eigen-problem, K v = w M v, using a matrix iteration approach
	with shifting.

	@param n number of degrees of freedom
	@param m number of required modes
*/
void stodola(
	float **K, float **M
	, int n, int m
	, float *w, float **V, float tol, float shift, int *iter, int *ok
);

#endif /* FRAME_EIG_H */

