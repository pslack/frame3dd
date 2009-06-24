
FRAME3DD - README for Windows users
-----------------------------------

This README file accompanies the binary packaged version of FRAME3DD for
the Windows platform.

See 'README.txt' for general information.

NOTE: FRAME3DD is a *command line* program. It processes an input text file
and adds output data at the end of that file.
It is best to run FRAME3DD from within a Command Prompt window.

NOTE: FRAME3DD input files, such as those in the 'examples' directory
currently include places where the file-paths for certain output files
will be written. In some cases, the locations will not exist on Windows
systems, so you will need to edit these input files before running the
program.

In order to keep your workspace 'clean', it is recommended to run
FRAME3DD on windows as follows:

  mkdir %HOMEPATH%\frame3dd
  cd %PROGRAMFILES%\FRAME3DD
  xcopy "examples\*" "%HOMEPATH%\frame3dd"
  cd %HOMEPATH%\frame3dd
  set PATH=%PATH%;%PROGRAMFILES%\FRAME3DD
  set FRAME3DD_OUTDIR=%HOMEPATH%\frame3dd
  frame3dd exA.3dd

After running FRAME3DD in this way, all your output files will be located 
in the the directory

c:\Documents and Setttings\yourname\frame3dd

(providing you have edited the input files so that they do not specify
absolute paths).

NOTES ON THE VIEWER PROGRAM
---------------------------

This installer includes a viewer for Microstran .ARC files. For this
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
