Name:		frame3dd
Summary:	Structural analysis of 2D/3D frames

# NOTE: IF THIS FILE IS NAMED 'frame3dd.spec' THEN DO NOT EDIT IT AS
# YOUR CHANGES WILL BE OVERWRITTEN. The 'frame3dd.spec' file is auto-generated
# from the 'frame3dd.spec.in' file; you should make your changes to this latter
# file instead.

# This version number is filled in automatically when you run 'scons dist'.
# You should update it in the 'SConstruct' file, rather than here.
Version:	0.20090417

# Use release 0.* so that other users can do patch releases with a higher number
# and still have the update occur automatically.
Release:	0%{?dist}

License:	GPL
Group:		Applications/Engineering
Source:		%{name}-%{version}.tar.bz2
URL:		http://frame3dd.sourceforge.net/index.html

Prefix:		%{_prefix}
Packager:	John Pye
Vendor:		Henri Gavin
BuildRoot:	%{_tmppath}/%{name}-%{version}-root

BuildRequires: gcc-c++
BuildRequires: scons >= 0.96.92
BuildRequires: SoQt-devel, Coin2-devel

%description
Free software for static and dynamic structural analysis of 2D and 3D frames 
and trusses with elastic and geometric stiffness. Computes the static
deflections, reactions, internal element forces, natural frequencies, mode
shapes and modal participation factors of two- and three- dimensional elastic
structures using direct stiffness and mass assembly.  Graphical output and
mode shape animation via Gnuplot version 4.0.

Requires: gnuplot, SoQT, Coin2.

%prep
%setup -q

%build
rm -rf %{buildroot}
scons %{?_smp_mflags} CXX="%{?ccache} g++" CC="%{?ccache} gcc" \
	INSTALL_ROOT=%{buildroot} \
	INSTALL_PREFIX=%{_prefix}

%install
scons %{?_smp_mflags} install \
	INSTALL_ROOT=%{buildroot} \
	INSTALL_PREFIX=%{_prefix}

%clean
rm -rf %{buildroot} 

%files
%defattr(-, root, root)
%doc doc/*
%{_bindir}
%{_libdir}/lib*.so*
%{_datadir}/%{name}

%changelog
#
# ChangeLog is now maintained in ChangeLog.txt. Make your
# changes there.
#
* Fri Apr 17 2009 John Pye <john.pye@anu.edu.au> 0.20090417
- Removed LinuxCommand PDF link in documentation.
- Removed little cmd prompt icons from documentation.
- New release on SF.net before new work on data structures commences.
- Reformatted changelog to PC format (may cause probs with RPM?)
- Removed wordy usage output from changelog (doesn't fit RPM formatting
requirements).

* Tue Mar 31 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090331
- Fixed bug in the Matlab interface function frame_3dd.m ... 
...no long printing the depricated anlyz variable value
- Frame3DD now writes the stiffness matrix to a file named "Ks" 
for each analysis. Line 460 of main.c ... after Newton-Raphson iterations 
for geometric nonlinear analysis, if such an analysis is to be performed.  
...for compatability with Matlab interface function frame_3dd.m
- re-ran example files

* Thu Mar 5 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090305 
- Remove "anlyz" from input data file format because the "-c"
command line option is now an easier and better way to specify
"data check only". Updated code, examples and documentation.
- Fixed bug related to "-q" flag and verbose output on line 1274 of frame3dd_io.c
- Added checks related to incorrect command-line arguments.
- Added an evaluation/interpretation of RMS relative equilibrium precision
  within the code and updated the documentation.

* Wed Mar 4 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090304 
- Fixed fprintf format character in save_ivector()
- Re-ran examples

* Wed Mar 4 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090304 
- Fixed bug in multi-load case nonlinear analysis in main.c line 410
  ... Iteration termination criteria must be reset at the start of each iteration.
- Change snprintf to sprintf in frame3dd_io.c for DJGPP compatability
- Added #include <time.h> in nrutil.c for DJGPP compatability

* Wed Mar 4 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090304 
- Added matrix condensation option to the command line
- Improved command line interface (now "frame3dd infile outfile" is ok)
- Improved user-manual.html regarding command-line options with examples
- Command-line help is (trimmed)

* Wed Mar 4 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090304
- Implement command line parsing using the getopt function
  ... using getopt for ease of portability.
- Added functions to frame3dd_io.c:
  parse_options() --- calls getopt
  display_help() 
  display_usage()
  display_version()
- Updated documentation with command-line syntax ... doc/user-manual.html

* Tue Mar 3 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090303
- Implement command line parsing using argtable2 package

* Mon Mar 2 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090302
- Usage change from "frame3dd InputData.frm" to "frame3dd
  InputData.frm OutputData.out"
- output information regarding number of loading types is now more clear
- updated examples B and E with trapezoidal loads
- updated documentation and README
- updated Matlab frame_3dd.m

* Thu Feb 27 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090227
- Implement trapezoidally-distributed loads over partial distances along frame elements

* Tue Feb 10 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090210
- Input and Output are now in separate files.  The Output filename is automatically generated from the Input file.  The Output file recapitualtes the input data before writing the output data.   
- Removed plot file, mesh file, and mode file names from input data file These file names are now automatically generated.  The path to these file names is also automatically generated according to the OS (Win32 or Linux/Unix/OSX etc). 
- Removed file names for the plot file, mesh file, and mode file from the examples. 
- Changed output_file_location to output_path
- Moved calls to output_path (for the /tmp or "TEMP" path) from main.c to frame3dd_io.c
- Changed  tpath to temppath
- Changed frame3dd.cln file name to frame3dd.frm ... clean file name
- All instances of a double quote character are ignored in reading input files
- Updated documentation.

* Mon Feb 09 2009 John Pye <john.pye@anu.edu.au> 0.20090209
- Added support for FRAME3DD_OUTDIR as location of output files
- Updated documentation
- Added datestamp to documentation (a bit kludgey)
- More comments on use of temp directory and output files.
- Changed example files to remove absolute /tmp paths.
- Changelog converted to separate file.
- Added missing documentation images to Windows installer.

* Mon Feb 02 2009 John Pye <john.pye@anu.edu.au> 0.20090202
- Fixed problem with use of /tmp on Windows
- Upload binaries to SF.net
- Add comments to README-win32.txt

* Thu Jan 29 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090129
- Completed migration of website to frame3dd.sourceforge.net

* Sun Jan 4 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090104
- Input data order reorganized: joints, reactions, members, loads
- Updated TODO.txt, website, examples, and documentation. 

* Thu Jan 1 2009 Henri Gavin <henri.gavin@duke.edu> 0.20090101
- Implemented .CSV spreadsheet support
- Eliminated some unneccessary vectors  "_lc"
- Updated TODO.txt, website, examples, and documentation. 

* Wed Dec 31 2008 Henri Gavin <henri.gavin@duke.edu> 0.20081231
- Moved rel_norm() from ldl_dcmp.c to frame3dd.c
- User input variables are now float instead of double
  ... except for joint location variable (xyz) which is still type vec3

* Tue Dec 30 2008 Henri Gavin <henri.gavin@duke.edu> 0.20081230
- The Input Data file format changed in two ways as follows:
  > Reaction data are listed just after the Member data.
  > Prescribed displacements are listed separately for each load case.
- Added comments to src/*.h files
- Changed "nM" variable name to "nB" ... number of beam elements
- Changed "modes" variable name to "nM" ... number of dynamic modes
- Changed "pos" variable name to "xyz" ... joint locations
- Removed many #include "abc.h" lines from .c source files
  without affecting the ability to execute or to compile warning-free via  ... 
    gcc -O -Wall -o frame3dd main.c frame3dd.c frame3dd_io.c ldl_dcmp.c lu_dcmp.c coordtrans.c eig.c nrutil.c -lm
- Updated TODO.txt, website, examples, and documentation.

* Fri Dec 12 2008 Henri Gavin <henri.gavin@duke.edu> 0.20081212
- The file frame3dd.cln is now written to /tmp/

* Thu Dec 11 2008 Henri Gavin <henri.gavin@duke.edu> 0.20081211
- Fixed bugs in writing Matlab m-file output.
- Documentation updated.
- The name of the outout m-file will not conflict with other file names.
- Renamed itoa function to my_itoa for portability reasons.

* Mon Dec 08 2008 Henri Gavin <henri.gavin@duke.edu> 0.20081208
- An error message is now displayed if nL < 1
- Updated documentation.

* Mon Dec 01 2008 Henri Gavin <henri.gavin@duke.edu> 0.20081201
- Removed the redefine of float to double, now all floating point
  variables are defined as doubles.
  In a future version variables that can actually be floats, such
  as variables read from input data files, will be re-defined as floats.
- Support for vec3 is retained througout.
- Renamed frm_io.c to frame3dd_io.c
- Renamed frm_io.h to frame3dd_io.h
- Updated Sconscript with these file name changes
- Changed getline_no_comment in frame3dd_io.c to support
  data files with comma-delimited values , with anticipated support
  for .CSV files

* Mon Oct 13 2008 John Pye <john@curioussymbols.com> 20081013
- Incorporated Henri's recent changes back into SF.net repository
- Changed URLs to http://frame3dd.sf.net (which currently redirects to
  Henri's personal homepage, but proposed to move to SF.net in future)
- Merged support for 'vec3' type in FRAME3DD, first step towards refactoring
  to some simpler API definitions.

* Tue Sep  9 2008 Henri Gavin <henri.gavin@duke.edu> 20080909
- Added multiple load case capability
- Updated web site http://www.duke.edu/~hpgavin/frame/
  with revised examples and revised instructions for the
  multiple load case capability

* Mon Mar 17 2008 John Pye <john@curioussymbols.com>
- Renamed to 'frame3dd' inline with new SF.net project name.

* Fri Mar 14 2008 Henri Gavin <henri.gavin@duke.edu> 20080314
- Changed name of project from FRAME to FRAME3DD
- Modified content of project web site ...
   http://www.duke.edu/~hpgavin/frame/
  and doc/user-manual.html to reflect change in name. The URL has not changed.
- Updated source files to reflect name change 
- Renamed src/frame.h to src/frame3dd.h
- Renamed src/frame.c to src/frame3dd.c
- Renamed frame.h to frame3dd.h
- Renamed frame.spec to frame3dd.spec
- Changed License from GPL version 2 to GPL version 3
- Changed license text on FRAME source code to GPL version 3
- Updated source files to reflect license change 
- Created Sourceforge project under the name FRAME3DD
- Files in src/viewer and src/microstran left unchanged.

* Mon Jan 15 2008 John Pye <john@curioussymbols.com> 20071220
- Updated to new version received from H Gavin by email.

* Wed Jun 28 2007 John Pye <john@curioussymbols.com> 20070301
- Start of Frame3DD SourceForge project

* Fri Jan 01 1993 Henri Gavin <henri.gavin@duke.edu> 19930101
- initiation of program at the University of Michigan



