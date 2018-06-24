export PARENT_URL=http://simvascular.stanford.edu/downloads/public/simvascular/externals/src/originals/

mkdir Originals
pushd Originals

#mkdir -p tcltk
#pushd tcltk
#wget $PARENT_URL/tcltk/tcl8.5.18-src.tar.gz
#wget $PARENT_URL/tcltk/tcl8.6.4-src.tar.gz
#wget $PARENT_URL/tcltk/tcllib-1.17.tar.gz
#wget $PARENT_URL/tcltk/tk8.5.18-src.tar.gz
#wget $PARENT_URL/tcltk/tk8.6.4-src.tar.gz
#wget $PARENT_URL/tcltk/tklib-0.6.tar.tgz
#popd

mkdir -p vtk
pushd vtk
wget $PARENT_URL/vtk/VTK-6.2.0.tar.gz
#wget $PARENT_URL/vtk/VTK-6.3.0.tar.gz
popd

popd
