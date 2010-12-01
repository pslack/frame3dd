/*       
 * This file is part of FRAME3DD: 
 * Static and dynamic structural analysis of 2D & 3D frames and trusses
 * with elastic and geometric stiffness.
 * ---------------------------------------------------------------------------
 * http://frame3dd.sourceforge.net/
 * ---------------------------------------------------------------------------
 * Copyright (C) 1992-2010  Henri P. Gavin
 *
 * FRAME3DD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FRAME3DD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FRAME3DD.  If not, see <http://www.gnu.org/licenses/>.
*//** @file
*/

#include <stdio.h>
#include <stdlib.h>

#define ANSI_SYS	1	/*  compile for ANSI_SYS driver; 0: don't */
// ... requires ANSI.SYS and the line   DEVICE = C:\ANSI.SYS  in  C:\CONFIG.SYS
// #define MAXL 128


/* ---------------------------------------------------------------------------
COLOR - change color on the screen ... 
 Screen   Color  Scheme  : 0 = white on black, 1 = bright
 first digit= 3  for text color          first digit= 4  for  background color
 second digit codes:     1=red, 2=green, 3=gold, 4=blue, 5=purple, 6=lght blue
---------------------------------------------------------------------------- */
void color ( const int colorCode );	/*  change the screen color      */


/* ---------------------------------------------------------------------------
TEXTCOLOR - change color of text and background
 tColor : text color : one of 'k' 'r' 'g' 'y' 'b' 'm' 'c' 'w'
 bColor : back color : one of 'k' 'r' 'g' 'y' 'b' 'm' 'c' 'w'
 nbf    : 'n' = normal, 'b' = bright/bold, 'f' = faint
 uline  : 'u' = underline
 http://en.wikipedia.org/wiki/ANSI_escape_code
--------------------------------------------------------------------------- */
void textColor ( const char tColor, const char bColor, const char nbf, const char uline );


/* ---------------------------------------------------------------------------
ERRORMSG -  write a diagnostic error message in color
---------------------------------------------------------------------------- */
void errorMsg ( const char *errString );


/*  -------------------------------------------------------------------------
OPENFILE  -  open a file or print a diagnostic error message 
---------------------------------------------------------------------------- */
FILE *openFile (const char *path, const char *fileName, const char *mode );


/* ---------------------------------------------------------------------------
SCANLINE -  scan through a line until a 'a' is reached, like getline() 3feb94
---------------------------------------------------------------------------- */
int scanLine ( FILE *fp, int lim, char *s, const char a );


/* ---------------------------------------------------------------------------
SCANLABEL -  scan through a line until a '"' is reached, like getline()
---------------------------------------------------------------------------- */
int scanLabel ( FILE *fp, int lim, char *s, const char a );


