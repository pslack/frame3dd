# !bin/bash
#
# dist_src.sh
# assemble .ZIP files for Frame3dd source code distribution release

export VERSION=20091022
echo $VERSION

# clean out prior distribution files
echo "cleaning out prior distribution files ... "
rm -rf dist/Frame3DD

# make directory structure for distribution .ZIP files in the ./dist directory
echo "creating directory structure for new distribution ... "
mkdir dist/Frame3DD
mkdir dist/Frame3DD/src
mkdir dist/Frame3DD/src/microstran

# copy documentation
echo "copying source code ... "
cp --preserve=mode,timestamps src/*.c               dist/Frame3DD/src/.
cp --preserve=mode,timestamps src/*.h               dist/Frame3DD/src/.
cp --preserve=mode,timestamps src/microstran/vec3.h dist/Frame3DD/src/microstran/.
cp --preserve=mode,timestamps src/microstran/config.h dist/Frame3DD/src/microstran/.

# copy ChangeLog, license, and README files
cp --preserve=mode,timestamps ChangeLog.txt         dist/Frame3DD/.
cp --preserve=mode,timestamps LICENSE.txt           dist/Frame3DD/.
cp --preserve=mode,timestamps README.txt            dist/Frame3DD/.
cp --preserve=mode,timestamps README-win32.txt      dist/Frame3DD/.

# make a Frame3DD/temp directory
mkdir dist/Frame3DD/temp

# assemble the .zip file
echo "assembling .zip file ... "
cd dist
# zip -r Frame3DD_$(date +%Y%m%d)_src.zip Frame3DD/*
zip -r Frame3DD_$(echo $VERSION)_src.zip Frame3DD/*

echo "Frame3DD .zip archive complete"

rm -rf Frame3DD

# ----------------------------------------------------------- dist_src.sh
# Henri P. Gavin 2009.10.20
# updated 2009.10.25, 2009.10.27
