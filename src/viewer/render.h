/*	FRAME3DD: Static and dynamic structural analysis of 2D & 3D frames and trusses
	Copyright (C) 2007-2008 John Pye

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
*/
#ifndef FRAME_RENDER_H
#define FRAME_RENDER_H

#include <microstran/sections.h>

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/SbVec3f.h>

#include <vector>

extern const SbColor WHITE;
extern const SbColor RED;
extern const SbColor GREEN;
extern const SbColor YELLOW;
extern const SbColor BLUE;
extern const SbColor ORANGE;
extern const SbColor PURPLE;
extern const SbColor MAGENTA;
extern const SbColor CYAN;
extern const SbColor LIME;

SoSeparator *text(const SbVec3f &left, const char *str, const SbColor &c=WHITE);
SoSeparator *sphere(const SbVec3f &C, const double &r=0.1, const SbColor &c=YELLOW);
SoSeparator *cylinder(const SbVec3f &A, const SbVec3f &B, const double &radius=0.1, const SbColor &c=YELLOW);
SoSeparator *cone(const SbVec3f &A, const SbVec3f &B, const double &r, const SbColor &c = YELLOW);
SoSeparator *axes(const double &size=1.0, double thickness=0.0, bool labelled=true);
SoSeparator *arrow(const SbVec3f &A, const SbVec3f &B, const SbColor &c=YELLOW, const char *label=NULL, double thickness=0.05);
SoSeparator *face(const SbVec3f &n, const std::vector<SbVec3f> &vertices, const SbColor &c=YELLOW);

SoSeparator *prism(const SbVec3f &A, const SbVec3f &B, const section_outline &outline, const SbColor &c, const SbVec3f &O);

#endif

