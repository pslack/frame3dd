# !bin/bash
#
# zipdist.sh
# assemble .ZIP files for Frame3dd executable and source distribution release

export VERSION=20100105
echo $VERSION

# clean out prior distribution files
echo "cleaning out prior distribution files ... "
rm -rf dist/Frame3DD

# make directory structure for distribution .ZIP files in the ./dist directory
echo "creating directory structure for new distribution ... "
mkdir dist/Frame3DD
mkdir dist/Frame3DD/doc
mkdir dist/Frame3DD/doc/img
mkdir dist/Frame3DD/examples
mkdir dist/Frame3DD/matlab

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

# copy example files
echo "copying example files ... "
cp --preserve=mode,timestamps examples/ex*.3dd      dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/ex*.out      dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/saveplot     dist/Frame3DD/examples/.
cp --preserve=mode,timestamps examples/saveplot_w32 dist/Frame3DD/examples/.

# copy matlab files
echo "copying matlab files ... "
cp --preserve=mode,timestamps matlab/*              dist/Frame3DD/matlab/.

# make an empty temp directory
echo "making an empty Frame3DD/temp directory ... "
mkdir dist/Frame3DD/temp/

# assemble the .zip file
echo "assembling .zip file ... "
cd dist                                           # change to trunk/dist
#zip -r Frame3DD_$(date +%Y%m%d).zip Frame3DD/*
zip -r Frame3DD_$(echo $VERSION).zip Frame3DD/*

# make copies of .zip files for each operating system
#cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_linux.zip
#cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_osx34.zip
#cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_osx56.zip
#cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_win32.zip

cp Frame3DD_$(echo $VERSION).zip Frame3DD_$(echo $VERSION)_DukeOIT.zip
cp Frame3DD_$(echo $VERSION).zip Frame3DD_$(echo $VERSION)_linux.zip
cp Frame3DD_$(echo $VERSION).zip Frame3DD_$(echo $VERSION)_osx34.zip
cp Frame3DD_$(echo $VERSION).zip Frame3DD_$(echo $VERSION)_osx56.zip
cp Frame3DD_$(echo $VERSION).zip Frame3DD_$(echo $VERSION)_win32.zip
cp Frame3DD_$(echo $VERSION).zip Frame3DD_$(echo $VERSION)_src.zip

#rm Frame3DD_$(date +%Y%m%d).zip
rm Frame3DD_$(echo $VERSION).zip

# add Linux executable
echo "adding Linux executable ... "
cp --preserve=mode,timestamps ../build/frame3dd        Frame3DD/.
# zip Frame3DD_$(date +%Y%m%d)_linux.zip Frame3DD/frame3dd
zip Frame3DD_$(echo $VERSION)_linux.zip Frame3DD/frame3dd
rm Frame3DD/frame3dd

# add Duke OIT executable  
echo "adding OS X executable ... "
cp --preserve=mode,timestamps ../build/frame3dd_oit     Frame3DD/frame3dd
zip Frame3DD_$(echo $VERSION)_DukeOIT.zip Frame3DD/frame3dd
rm Frame3DD/frame3dd

# add OS X 10.3 10.4 10.5 10.6 executables  
echo "adding OS X executable ... "
cp --preserve=mode,timestamps ../build/frame3ddosx34     Frame3DD/frame3dd
zip Frame3DD_$(echo $VERSION)_osx34.zip Frame3DD/frame3dd
cp --preserve=mode,timestamps ../build/frame3ddosx56       Frame3DD/frame3dd
zip Frame3DD_$(echo $VERSION)_osx56.zip Frame3DD/frame3dd
rm Frame3DD/frame3dd

# add Windows executable
echo "adding Windows executable ... "
cp --preserve=mode,timestamps ../build/frame3dd.exe    Frame3DD/.
#zip Frame3DD_$(date +%Y%m%d)_win32.zip Frame3DD/frame3dd.exe
zip Frame3DD_$(echo $VERSION)_win32.zip Frame3DD/frame3dd.exe
rm Frame3DD/frame3dd.exe

echo "Frame3DD executable .zip archive complete"

cd ..

# add Frame3dd source code to distribution release

# make directory structure for distribution .ZIP files in the ./dist directory
echo "creating source code directory structure for new distribution ... "
mkdir dist/Frame3DD/src
mkdir dist/Frame3DD/src/microstran

# copy source code
echo "copying source code ... "
cp --preserve=mode,timestamps src/*.c               dist/Frame3DD/src/.
cp --preserve=mode,timestamps src/*.h               dist/Frame3DD/src/.
cp --preserve=mode,timestamps src/Makefile          dist/Frame3DD/src/.
cp --preserve=mode,timestamps src/microstran/vec3.h dist/Frame3DD/src/microstran/.
cp --preserve=mode,timestamps src/microstran/config.h dist/Frame3DD/src/microstran/.

# assemble the .zip file
echo "adding source code to Frame3DD_*_src.zip file ... "
cd dist
# zip -r Frame3DD_$(date +%Y%m%d)_src.zip Frame3DD/src/*
zip -r Frame3DD_$(echo $VERSION)_src.zip Frame3DD/src/*

echo "Frame3DD source code .zip archive complete"

rm -rf Frame3DD

# To upload distribution files to sourceforge.net ...
#
# rsync -uav *.zip hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20100105/.
#
# rsync -uav *.bz2  hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20100105/.

# ----------------------------------------------------------- zipdist.sh
# Henri P. Gavin  2009.10.20
# updated: 2009.10.22, 2009.10.25, 2009.10.27, 2009.10.29, 2010.01.05, 2010.12.1
