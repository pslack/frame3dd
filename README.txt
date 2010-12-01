FRAME3DD: a frame analysis program
-------------------------------

FRAME3DD is a program for the static and dynamic structural analysis of
two- and three-dimensional frames and trusses with elastic and 
geometric stiffness.

FRAME3DD reads an Input Data file, containing joint coordinates, member
geometry, material moduli, restrained joints, prescribed displacements,
load information, and optionally, mass information if a modal analysis
is to be carried out.

FRAME3DD appends results to an Output Data File, The Output Data 
recapitulates the input information, gives joint displacements in global 
coordinates, member end-forces in local coordinates, reactions in 
global coordinates, and natural frequencies and mode shapes in global 
coordinates.

For usage information, see the User Manual, included with both the source
code and binary installer packages.

Additional output files from FRAME3DD are directed by default to the temporary
file folder (according to the %TEMP% env var on Windows, or /tmp on other
platforms). You can specify an alternative location for output files by setting
the environment variable FRAME3DD_OUTDIR.

Note that you require write permissions for the input file in order that
FRAME3DD can append its calculation results to that file.

---

FRAME3DD is free software; you may redistribute it and/or modify it under 
the terms of the GNU General Public License (GPL) as published by the
Free Software Foundation. FRAME3DD is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License (GPL) for details. 

---

You should have obtained this distribution from SourceForge.net 

svn co https://frame3dd.svn.sourceforge.net/svnroot/frame3dd frame3dd 

---

Included in this distribution:

LICENSE.txt  ... GNU GPL license
README.txt   ... this file

SConstruct   ... optionally used for compiling and building distribution
                 using scons ... http://www.scons.org/
frame3dd.spec... build specification for scons, includes changelog
save.plt     ... Gnuplot macro for saving Gnuplot plots as PostScript

doc/         ... documentation (in HTML format)
examples/    ... examples of FRAME input/output files
scons/       ... used by scons to build distribution
src/         ... C source code 

to compile and build without scons ...

gcc -O -Wall -o frame3dd main.c frame3dd.c frame3dd_io.c ldl_dcmp.c lu_dcmp.c coordtrans.c eig.c hpgUtils.c nrutil.c -lm

---

FRAME3DD author

(c) 1992-2008 Henri Gavin - Associate Professor 
Department of Civil and Environmental Engineering
Edmund T. Pratt School of Engineering
Duke University - Box 90287, Durham, NC 27708-0287

mailto:Henri.Gavin@Duke.edu - tel: 919-660-5201 - fax: 919-660-5219 

---

Packaging of FRAME3DD for distribution on Linux and Windows,
as well as implementation of the Microstran file parsers
is by John Pye, 2007-2009.

John Pye
Department of Engineering
Australian National University
mailto:john.pye@anu.edu.au
http://pye.dyndns.org/

---

