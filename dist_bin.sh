# !bin/bash
#
# dist_zip.sh
# assemble .ZIP files for Frame3dd binary distribution release

# clean out prior distribution files
echo "cleaning out prior distribution files ... "
rm -rf dist/Frame3DD

# make directory structure for distribution .ZIP files in the ./dist directory
echo "creating directory structure for new distribution ... "
mkdir dist/Frame3DD
mkdir dist/Frame3DD/doc
mkdir dist/Frame3DD/doc/img
mkdir dist/Frame3DD/doc/Matlab
mkdir dist/Frame3DD/examples

# copy documentation
echo "copying documentation ... "
cp --preserve=mode,timestamps doc/*.html            dist/Frame3DD/doc/.
cp --preserve=mode,timestamps doc/template.3dd      dist/Frame3DD/doc/.
cp --preserve=mode,timestamps doc/img/*.png         dist/Frame3DD/doc/img/.
cp --preserve=mode,timestamps doc/img/*.jpg         dist/Frame3DD/doc/img/.
cp --preserve=mode,timestamps doc/img/*.gif         dist/Frame3DD/doc/img/.
cp --preserve=mode,timestamps doc/Matlab/*          dist/Frame3DD/doc/Matlab/.
cp --preserve=mode,timestamps doc/LinuxTerminal.pdf dist/Frame3DD/doc/.	

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

# copy executables
echo "copying executables ... "
cp --preserve=mode,timestamps build/frame3dd.exe    dist/Frame3DD/.
cp --preserve=mode,timestamps build/frame3dd        dist/Frame3DD/.
cp --preserve=mode,timestamps build/frame3ddosx     dist/Frame3DD/.

# assemble the .zip file
echo "assembling .zip file ... "
cd dist
zip -r Frame3DD_$(date +%Y%m%d).zip Frame3DD/*

# delete non-os executables from each .zip distribution
cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_linux.zip
cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_osx.zip
cp Frame3DD_$(date +%Y%m%d).zip Frame3DD_$(date +%Y%m%d)_win32.zip

zip -d Frame3DD_$(date +%Y%m%d)_linux.zip Frame3DD/frame3dd.exe
zip -d Frame3DD_$(date +%Y%m%d)_linux.zip Frame3DD/frame3ddosx

zip -d Frame3DD_$(date +%Y%m%d)_osx.zip Frame3DD/frame3dd.exe
zip -d Frame3DD_$(date +%Y%m%d)_osx.zip Frame3DD/frame3dd

zip -d Frame3DD_$(date +%Y%m%d)_win32.zip Frame3DD/frame3ddosx
zip -d Frame3DD_$(date +%Y%m%d)_win32.zip Frame3DD/frame3dd

rm Frame3DD_$(date +%Y%m%d).zip

echo "Frame3DD .zip archive complete"

# Uploading distribution files to sourceforge.net ...
#
# scp Frame3DD_20091020_linux.zip hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20091020/.
#
# scp Frame3DD_20091020_win32.zip hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20091020/.
#
# scp Frame3DD_20091020_osx.zip hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20091020/.
#
# scp Frame3DD_20091020_src.zip hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20091020/.
#
# scp frame3dd-0.20091020.tar.bz2  hpgavin,frame3dd@frs.sourceforge.net:/home/frs/project/f/fr/frame3dd/frame3dd/0.20091020/.

# ----------------------------------------------------------- dist_bin.sh
# Henri P. Gavin  2009.10.20
# updated: 2009.10.22
