#
# unpack all of the source code
#

source Scripts/untar-unzip-source-all.sh
mkdir -p tmp

#
# make build scripts
#


#  tcl/tk 8.6
#sed -f CompileScripts/sed-script-x64_linux-options-gcc.sh CompileScripts/tcl-linux-generic.sh > tmp/compile.make.tcl.gcc.sh
#chmod a+rx ./tmp/compile.make.tcl.gcc.sh

# vtk
sed -f CompileScripts/sed-script-x64_linux-options-gcc.sh CompileScripts/compile-cmake-vtk-generic.sh > tmp/compile.cmake.vtk.gcc.sh
chmod a+rx ./tmp/compile.cmake.vtk.gcc.sh

# create script to create tar files
sed -f CompileScripts/sed-script-x64_linux-options-gcc.sh Scripts/create-archives-generic.sh > tmp/create-archives-all.gcc.sh
chmod a+rx ./tmp/create-archives-all.gcc.sh

# create script to create zip files
sed -f CompileScripts/sed-script-x64_linux-options-gcc.sh Scripts/tar-to-zip-all.sh > tmp/tar-to-zip-all.gcc.sh
chmod a+rx ./tmp/tar-to-zip-all.gcc.sh

#
# compile code
#

#  tcl/tk 8.6
#./tmp/compile.make.tcl.gcc.sh >& ./tmp/stdout.tcl.txt

# vtk
./tmp/compile.cmake.vtk.gcc.sh >& ./tmp/stdout.vtk.gcc.txt

#
# create tar files for distrution
#

./tmp/create-archives-all.gcc.sh >& ./tmp/stdout.create-archives-all.gcc.txt

./tmp/tar-to-zip-all.gcc.sh >& ./tmp/stdout.tar-to-zip-all.gcc.txt
