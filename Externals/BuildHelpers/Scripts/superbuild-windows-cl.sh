#
# unpack all of the source code
#

source Scripts/untar-unzip-source-all.sh
mkdir -p tmp

#
# make build scripts
#

#  tcl/tk 8.6
#sed -f CompileScripts/sed-script-x64_cygwin-options-cl.sh CompileScripts/tcl-windows-generic.sh > tmp/compile.make.tcl.cl.sh
#chmod a+rx ./tmp/compile.make.tcl.cl.sh

# vtk
sed -f CompileScripts/sed-script-x64_cygwin-options-cl.sh CompileScripts/compile-cmake-vtk-generic.sh > tmp/compile.cmake.vtk.cl.sh
chmod a+rx ./tmp/compile.cmake.vtk.cl.sh

# create script to create tar files
sed -f CompileScripts/sed-script-x64_cygwin-options-cl.sh Scripts/create-archives-generic.sh > tmp/create-archives-windows.cl.sh
chmod a+rx ./tmp/create-archives-windows.cl.sh

# create script to create zip files
sed -f CompileScripts/sed-script-x64_cygwin-options-cl.sh Scripts/tar-to-zip-all.sh > tmp/tar-to-zip-all.windows.cl.sh
chmod a+rx ./tmp/tar-to-zip-all.windows.cl.sh

#
# compile code
#

#  tcl/tk 8.6
#./tmp/compile.make.tcl.cl.sh >& ./tmp/stdout.tcl.txt

# vtk
./tmp/compile.cmake.vtk.cl.sh >& ./tmp/stdout.vtk.cl.txt

#
# create tar files for distrution
#

./tmp/create-archives-windows.cl.sh >& ./tmp/stdout.create-archives-windows.cl.txt

./tmp/tar-to-zip-all.windows.cl.sh >& ./tmp/stdout.tar-to-zip-all.windows.cl.txt
