This folder contains a parser for Microstran .arc files. This is a text-file
export format supported by the commercial programs Microstran, Space Gass and
Multiframe.

TODO need documentation of this file format, if possible.

The parser is implemented using a pretty low-level pure-C method described here:
http://www.math.chalmers.se/~koen/ParserComboC/parser-combo-c.html

At this point, the parser code has been tested but hasn't yet been connected
to FRAME to build up an actual model structure that can be calculated.

