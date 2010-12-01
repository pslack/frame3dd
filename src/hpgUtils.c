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
#include "hpgUtils.h"
#include "common.h"

#define DEBUG 0


/* -----------------------------------------------------------------------------
COLOR - change color on the screen ... 
 Screen   Color  Scheme  : 0 = white on black, 1 = bright
 first digit= 3  for text color          first digit= 4  for  background color
 second digit codes:     1=red, 2=green, 3=gold, 4=blue, 5=purple, 6=lght blue
 http://en.wikipedia.org/wiki/ANSI_escape_code
------------------------------------------------------------------------------*/
void color ( const int colorCode )	/*  change the screen color      */
{
#if ANSI_SYS
	fprintf (stderr, "\033[%02dm", colorCode );
	(void) fflush(stderr);
#endif
	return;
}


/* ---------------------------------------------------------------------------
TEXTCOLOR - change color of text and background
 tColor : text color : one of 'k' 'r' 'g' 'y' 'b' 'm' 'c' 'w'
 bColor : back color : one of 'k' 'r' 'g' 'y' 'b' 'm' 'c' 'w'
 nbf    : 'n' = normal, 'b' = bright/bold, 'f' = faint
 uline  : 'u' = underline
 http://en.wikipedia.org/wiki/ANSI_escape_code
---------------------------------------------------------------------------- */
void textColor ( const char tColor, const char bColor, const char nbf, const char uline )
{
#if ANSI_SYS
	fprintf (stderr, "\033[%02d",0);// Control Sequence Introducer & reset		
	// background colors
	if ( bColor == 'k' ) fprintf (stderr, ";%02d", 40 ); // black
	if ( bColor == 'r' ) fprintf (stderr, ";%02d", 41 ); // red
	if ( bColor == 'g' ) fprintf (stderr, ";%02d", 42 ); // green
	if ( bColor == 'y' ) fprintf (stderr, ";%02d", 43 ); // yellow
	if ( bColor == 'b' ) fprintf (stderr, ";%02d", 44 ); // blue
	if ( bColor == 'm' ) fprintf (stderr, ";%02d", 45 ); // magenta
	if ( bColor == 'c' ) fprintf (stderr, ";%02d", 46 ); // cyan
	if ( bColor == 'w' ) fprintf (stderr, ";%02d", 47 ); // white

	// text colors
	if ( tColor == 'k' ) fprintf (stderr, ";%02d", 30 ); // black
	if ( tColor == 'r' ) fprintf (stderr, ";%02d", 31 ); // red
	if ( tColor == 'g' ) fprintf (stderr, ";%02d", 32 ); // green
	if ( tColor == 'y' ) fprintf (stderr, ";%02d", 33 ); // yellow
	if ( tColor == 'b' ) fprintf (stderr, ";%02d", 34 ); // blue
	if ( tColor == 'm' ) fprintf (stderr, ";%02d", 35 ); // magenta
	if ( tColor == 'c' ) fprintf (stderr, ";%02d", 36 ); // cyan
	if ( tColor == 'w' ) fprintf (stderr, ";%02d", 37 ); // white

//	printf(" tColor = %c   bColor = %c   nbf = %c\n", tColor, bColor, nbf );
	if ( nbf    == 'b' ) fprintf (stderr, ";%02d",  1 ); // bright
	if ( nbf    == 'f' ) fprintf (stderr, ";%02d",  2 ); // faint

	if ( uline == 'u' )  fprintf (stderr, ";%02d", 4 );  // underline

	fprintf (stderr,"m");		// Select Graphic Rendition (SGR)

	(void) fflush(stderr);
#endif
	return;
}


/* ---------------------------------------------------------------------------
ERRORMSG -  write a diagnostic error message in color
-----------------------------------------------------------------------------*/
void errorMsg ( const char *errString )
{
#if ANSI_SYS
	color(1); color(41); color(37);
#endif
	fprintf(stderr,"  %s  ", errString );
#if ANSI_SYS
	fflush(stderr);
	color(0);
#endif
	fprintf(stderr,"\n\n");
	return;
}


/* ---------------------------------------------------------------------------
OPENFILE  -  open a file or print a diagnostic error message 
-----------------------------------------------------------------------------*/
FILE *openFile ( const char *path, const char *fileName, const char *mode )
{
	FILE	*fp;
	char	pathToFile[MAXL], errMsg[MAXL];

	if (mode == 0)	return 0;

	sprintf(pathToFile,"%s%s", path, fileName );
#if DEBUG
	printf(" openFile ... file name = %s\n", pathToFile);
#endif
	if ((fp=fopen(pathToFile,mode)) == NULL ) { // open file 
		switch (*mode) {
		   case 'w':
			sprintf(errMsg,"%s%s","cannot write to file: ", pathToFile );
			break;
		   case 'r':
			sprintf(errMsg,"%s%s","cannot read from file: ", pathToFile );
			break;
		   case 'a':
			sprintf(errMsg,"%s%s","cannot append to file: ", pathToFile );
			break;
		   default:
			sprintf(errMsg,"%s%s","cannot open file: ", pathToFile );
		}
		errorMsg ( errMsg );
		exit(1);
	} else {
#if DEBUG
	printf(" openFile ... fp = %x\n", fp);
#endif
	
		return fp;
	}
}


/* ---------------------------------------------------------------------------
SCANLINE -  scan through a line until a 'a' is reached, like getline() 3feb94
-----------------------------------------------------------------------------*/
int scanLine ( FILE *fp, int lim, char *s, const char a ) 
{
       	int     c=0,  i=-1;

	while (--lim > 0 && (c=getc(fp)) != EOF && c != a)  s[++i] = c;
	s[++i]='\0';
	return i ;
}


/* ---------------------------------------------------------------------------
SCANLABEL -  scan through a line until a '"' is reached, like getline()
-----------------------------------------------------------------------------*/
int scanLabel ( FILE *fp, int lim, char *s, const char a )
{
       	int     c=0,  i=-1;

	while (--lim > 0 && (c=getc(fp)) != EOF && c != a)
		;			// scan to first delimitter char
	while (--lim > 0 && (c=getc(fp)) != EOF && c != a) 
		s[++i] = c;		// read the label between delimitters
	s[++i]='\0';
	return i ;
}


#undef DEBUG
