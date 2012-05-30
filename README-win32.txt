
FRAME3DD - README FOR WINDOWS USERS AND DEVELOPERS
--------------------------------------------------

This README file accompanies the binary packaged version of FRAME3DD for
the Windows platform.

See 'README.txt' for general information.


NOTES FOR WINDOWS USERS
-----------------------

The user manual (doc/user-manual.html) includes instructions for installing
and running Frame3DD under Windows.    Following these instructions will
install Frame3DD to your Desktop.   

FRAME3DD is a *command line* program. It processes an input text file
and adds output data at the end of that file.

It is best to run FRAME3DD from within a Command Prompt (or Console) window.



NOTES REGARDING ANSI.SYS and COLORS 
-----------------------------------

Frame3DD displays text in color using ANSI.SYS escape sequences.
Current versions (Windows 7, VISTA, etc) of the Microsoft
"Command Prompt" do not render ANSI.SYS escape sequences correctly  ... 
"Console windows in Windows versions based on NT (Windows NT 4.0, Windows 2000,
Windows XP, Windows Server 2003, Windows Vista and Windows Server 2008)
do not support ANSI Escape sequences at all."
http://en.wikipedia.org/wiki/ANSI_escape_code#Windows_and_DOS

Instead of showing text in color, the Microsoft "Command Prompt" program   
displays strange text sequences.   

If these text sequences become annoying to look at, you have two options.


ANSI.SYS option #1:

Download an ANSI capable console for Windows ... ANSICON ... from
  ... http://adoxa.110mb.com/ansicon/ ...
Extract the .ZIP file and copy the folder "x86" into the Desktop\Frame3DD\ 
Make a shortcut from Desktop/Frame3DD/x86/ansicon to the Desktop.
Run ANSICON instead of the Microsoft "Command Prompt"


ANSI.SYS option #2:

Recompile Frame3DD without the ANSI.SYS escape sequences by changing
line 28 of hpgUtils.h to ...

#define ANSI_SYS        0 

... before compiling.    




NOTES ON THE MICROSTRAN .ARC VIEWER PROGRAM
-------------------------------------------

The Windows installer includes a viewer for Microstran .ARC files. For this
viewer to run correctly, you must do a bit of work.

1. Install QT 4.3.3.
2. Install Coin3D 2.5.0
3. Install SoQT 1.3

The above installers are accessible via the following web page:
http://ascendwiki.cheme.cmu.edu/Binary_installers_for_Coin3d_and_SoQt_on_MinGW

4. Get the mingwm10.dll DLL file from the following archive, and copy it into
your c:\Program Files\FRAME3DD folder.

http://sourceforge.net/project/downloading.php?group_id=2435&filename=mingwrt-3.15.2-mingw32-dll.tar.gz

5. Set your PATH to contain the DLL folders for the packages in steps 1 
   to 3. Depending on where you installed everything, this could be:
   
   c:\QT\4.3.3\bin;c:\Program Files\Coin3D-2.5.0\bin;c:\Program Files\SoQt-1.4.1
   
After performing these steps, you should be able to double-click on an ARC 
file and it should display in a viewer window.

-- 
John Pye
http://pye.dyndns.org/ (<-- contact details here)
Jun 2009.
