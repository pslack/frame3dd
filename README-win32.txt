
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

-- 
John Pye
Feb 2009.
