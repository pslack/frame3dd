#!/bin/bash
#
# zipdist.sh
# assemble .ZIP files for Frame3dd executable and source distribution release

# NOTE: this version number should preferably be synchronised with the version
# number in file 'SConstruct'. Note the alternative packaging instructions in 
# the file 'PACKAGING.txt' which builds DEB, RPM and EXE packages, but nothing
# for Mac users :-( 
export VERSION=20140514+
echo $VERSION

# clean out prior distribution files
echo "cleaning out prior distribution files ... "
rm -rf dist/Frame3DD

# make directory structure for distribution .ZIP file in the ./dist directory
echo "creating directory structure for new distribution ... "
mkdir dist/Frame3DD
mkdir dist/Frame3DD/doc
mkdir dist/Frame3DD/doc/img
mkdir dist/Frame3DD/examples
mkdir dist/Frame3DD/matlab
mkdir dist/Frame3DD/linux
mkdir dist/Frame3DD/osx
mkdir dist/Frame3DD/windows

# make an empty temp directory
echo "making an empty Frame3DD/temp directory ... "
mkdir dist/Frame3DD/temp/

echo "creating source code directory structure for new distribution ... "
mkdir dist/Frame3DD/src
mkdir dist/Frame3DD/src/microstran

# copy documentation
echo "copying documentation ... "
cp --preserve=mode,timestamps doc/*.html            dist/Frame3DD/doc/.
cp --preserve=mode,timestamps doc/template.3dd      dist/Frame3DD/doc/.
cp --preserve=mode,timestamps doc/template.csv      dist/Frame3DD/doc/.
cp --preserve=mode,timestamps doc/img/*.png         dist/Frame3DD/doc/img/.
cp --preserve=mode,timestamps doc/img/*.jpg         dist/Frame3DD/doc/img/.
cp --preserve=mode,timestamps doc/img/*.gif         dist/Frame3DD/doc/img/.
#cp --preserve=mode,timestamps doc/LinuxTerminal.pdf dist/Frame3DD/doc/.	

# copy ChangeLog, license, and README files
cp --preserve=mode,timestamps ChangeLog.txt         dist/Frame3DD/.
cp --preserve=mode,timestamps LICENSE.txt           dist/Frame3DD/.
cp --preserve=mode,timestamps README.txt            dist/Frame3DD/.
cp --preserve=mode,timestamps README-win32.txt      dist/Frame3DD/.

# copy example files (text and .csv formats)
echo "copying example files ... "
cp --preserve=mode,timestamps examples/ex*.3dd      dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/ex*.out      dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/ex*.csv      dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/ex*_out.CSV  dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/saveplot     dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/saveplot_w32 dist/Frame3DD/examples/.

# copy matlab files
echo "copying matlab files ... "
cp --preserve=mode,timestamps matlab/*              dist/Frame3DD/matlab/.

# copy source code
echo "copying source code ... "
cp --preserve=mode,timestamps src/*.c               dist/Frame3DD/src/.
cp --preserve=mode,timestamps src/*.h               dist/Frame3DD/src/.
cp --preserve=mode,timestamps src/Makefile          dist/Frame3DD/src/.
cp --preserve=mode,timestamps src/microstran/vec3.h dist/Frame3DD/src/microstran/.
cp --preserve=mode,timestamps src/microstran/config.h dist/Frame3DD/src/microstran/.


# copy binary executable application
echo "copying Linux executable ... "
cp --preserve=mode,timestamps build/linux/frame3dd        dist/Frame3DD/linux/frame3dd
echo "copying OS X executable ... "
cp --preserve=mode,timestamps build/osx/frame3dd          dist/Frame3DD/osx/frame3dd
echo "copying Windows executable ... "
cp --preserve=mode,timestamps build/windows/frame3dd.exe  dist/Frame3DD/windows/frame3dd.exe
cp --preserve=mode,timestamps build/DukeOIT/frame3dd      dist/Frame3DD/DukeOIT/frame3dd

# assemble the .zip file
echo "assembling .zip file ... "
cd dist                                           # change to trunk/dist
zip -r Frame3DD_$(echo $VERSION).zip Frame3DD/*
#zip -r Frame3DD_$(date +%Y%m%d).zip Frame3DD/*

#rm Frame3DD_$(date +%Y%m%d).zip

echo "Frame3DD .zip archive complete"

rm -rf Frame3DD

# To upload distribution files to sourceforge.net ...
#
# rsync -uav Frame3DD*.zip hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/.
#
##rsync -uav *.bz2  hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20140514+/.

# ----------------------------------------------------------------- zipdist.sh
# Henri P. Gavin  2009-10-20
# updated: 2009-10-22, 2009-10-25, 2009-10-27, 2009-10-29, 2010-01-05, 2010-12-1, 2013-03-18, 2014-01-21, 2014-05-14, 2014-05-17, 2015-06-04
