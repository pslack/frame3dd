FRAME: a frame analysis program
-------------------------------

FRAME is a program for the static and dynamic structural analysis of
two- and three-dimensional frames and trusses with elastic and 
geometric stiffness.

FRAME reads an Input Data file, containing joint coordinates, member
geometry, material moduli, restrained joints, prescribed displacements,
load information, and optionally, mass information if a modal analysis
is to be carried out.

FRAME appends the Input Data with Output Data, resulting in a single
Input/Output file. The Output Data recapitulates the input information,
gives joint displacements in global coordinates, member end-forces in
local coordinates, reactions in global coordinates, and natural
frequencies and mode shapes in global coordinates.

FRAME is free software; you may redistribute it and/or modify it under 
the terms of the GNU General Public License (GPL) as published by the
Free Software Foundation. FRAME is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License (GPL) for details. 

---

Included in this distribution:

LICENSE.txt  ... GNU GPL license
README.txt   ... this file
SConstruct   ... optionally used for compiling and building distribution
                 using scons ... http://www.scons.org/
frame.spec   ... build specification for scons, includes changelog
saveplot     ... Gnuplot macro for saving Gnuplot plots as PostScript

doc/         ... documentation 
examples/    ... examples of FRAME input/output files
scons/       ... used by scons to build distribution
src/         ... C source code 

---

FRAME author:

(c) 1992-2007 Henri Gavin - Associate Professor 
Department of Civil and Environmental Engineering
Edmund T. Pratt School of Engineering
Duke University - Box 90287, Durham, NC 27708-0287

Henri.Gavin@Duke.edu - tel: 919-660-5201 - fax: 919-660-5219 

---

Packaged tarball and build script/installer by John Pye, 2007.
Department of Engineering
Australian National University
http://pye.dyndns.org/


