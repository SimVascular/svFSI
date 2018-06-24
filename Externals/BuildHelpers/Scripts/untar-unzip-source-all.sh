#
#  untar tcl/tk
#

#tar xvf Originals/tcltk/tcl8.5.18-src.tar.gz
#tar xvf Originals/tcltk/tk8.5.18-src.tar.gz  
#tar xvf Originals/tcltk/tcl8.6.4-src.tar.gz
#tar xvf Originals/tcltk/tk8.6.4-src.tar.gz
#tar xvf Originals/tcltk/tcllib-1.17.tar.gz
#tar xvf Originals/tcltk/tklib-0.6.tar.tgz

#
# move and rename tcl/tk
#

#mv tcl8.5.18 ../tcl-8.5.18
#mv tk8.5.18 ../tk-8.5.18
#mv tcl8.6.4 ../tcl-8.6.4
#mv tk8.6.4 ../tk-8.6.4
#mv tcllib-1.17 ../tcllib-1.17
#mv tklib-0.6 ../tklib-0.6
#source Patches/patch-source-tcltk-8.5.sh

#
# vtk
#

tar xvf Originals/vtk/VTK-6.2.0.tar.gz
mv VTK-6.2.0 ../vtk-6.2.0
source Patches/patch-source-vtk-6.2.sh
#tar xvf Originals/vtk/VTK-6.3.0.tar.gz
#mv VTK-6.3.0 ../vtk-6.3.0
#source Patches/patch-source-vtk-6.3.sh



