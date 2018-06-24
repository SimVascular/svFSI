#
# unpack all of the source code
#

source Scripts/untar-unzip-source-all.sh
mkdir -p tmp

#
# make build scripts
#


#  tcl/tk 8.6
#sed -f CompileScripts/sed-script-x64_mac_osx-options-clang.sh CompileScripts/tcl-mac_osx-generic.sh > tmp/compile.make.tcl.clang.sh
#chmod a+rx ./tmp/compile.make.tcl.clang.sh

# vtk
sed -f CompileScripts/sed-script-x64_mac_osx-options-clang.sh CompileScripts/compile-cmake-vtk-generic.sh > tmp/compile.cmake.vtk.clang.sh
chmod a+rx ./tmp/compile.cmake.vtk.clang.sh

# create script to create tar files
sed -f CompileScripts/sed-script-x64_mac_osx-options-clang.sh Scripts/create-archives-generic.sh > tmp/create-archives-mac_osx.clang.sh
chmod a+rx ./tmp/create-archives-mac_osx.clang.sh

# create script to create zip files
sed -f CompileScripts/sed-script-x64_mac_osx-options-clang.sh Scripts/tar-to-zip-all.sh > tmp/tar-to-zip-all.clang.sh
chmod a+rx ./tmp/tar-to-zip-all.clang.sh

#
# compile code
#

#  tcl/tk 8.6
#./tmp/compile.make.tcl.clang.sh >& ./tmp/stdout.tcl.txt

# vtk
./tmp/compile.cmake.vtk.clang.sh >& ./tmp/stdout.vtk.clang.txt

#
# create tar files for distrution
#

./tmp/create-archives-mac_osx.clang.sh >& ./tmp/stdout.create-archives-mac_osx.clang.txt

./tmp/tar-to-zip-all.clang.sh >& ./tmp/stdout.tar-to-zip-all.clang.txt
