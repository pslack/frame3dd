Name:		frame3dd
Summary:	Structural analysis of 2D/3D frames

# NOTE: IF THIS FILE IS NAMED 'frame3dd.spec' THEN DO NOT EDIT IT AS
# YOUR CHANGES WILL BE OVERWRITTEN. The 'frame3dd.spec' file is auto-generated
# from the 'frame3dd.spec.in' file; you should make your changes to this latter
# file instead.

# This version number is filled in automatically when you run 'scons dist'.
# You should update it in the 'SConstruct' file, rather than here.
Version:	0.20080911

# Use release 0.* so that other users can do patch releases with a higher number
# and still have the update occur automatically.
Release:	0%{?dist}

License:	GPL
Group:		Applications/Engineering
Source:		%{name}-%{version}.tar.bz2
URL:		http://www.duke.edu/~hpgavin/frame/index.html

Prefix:		%{_prefix}
Packager:	John Pye
Vendor:		Henri Gavin
BuildRoot:	%{_tmppath}/%{name}-%{version}-root

BuildRequires: gcc-c++
BuildRequires: scons >= 0.96.92

%description
Free software for static and dynamic structural analysis of 2D and 3D frames 
and trusses with elastic and geometric stiffness. Computes the static
deflections, reactions, internal element forces, natural frequencies, mode
shapes and modal participation factors of two- and three- dimensional elastic
structures using direct stiffness and mass assembly.  Graphical output and
mode shape animation via Gnuplot version 4.0.

Requires: gnuplot

%prep
%setup -q

%build
rm -rf %{buildroot}
scons %{?_smp_mflags} CXX="%{?ccache} g++" \
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
%{_datadir}

%changelog

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
- Initial version

