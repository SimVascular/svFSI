#/bin/bash  -f

# Copyright (c) 2015 Open Source Medical Software Corporation.
# All Rights Reserved.

REPLACEME_SV_SPECIAL_COMPILER_SCRIPT

export CC=CL

rm -Rf REPLACEME_SV_TOP_BIN_DIR_TCL
mkdir -p REPLACEME_SV_TOP_BIN_DIR_TCL
chmod u+w,a+rx REPLACEME_SV_TOP_BIN_DIR_TCL


cd ../REPLACEME_SV_TCL_DIR/win
nmake -f makefile.vc MACHINE=AMD64 INSTALLDIR="REPLACEME_SV_TOP_BIN_DIR_TCL" OPTS=msvcrt,static,threads hose
nmake -f makefile.vc MACHINE=AMD64 INSTALLDIR="REPLACEME_SV_TOP_BIN_DIR_TCL" OPTS=msvcrt,static,threads release
nmake -f makefile.vc MACHINE=AMD64 INSTALLDIR="REPLACEME_SV_TOP_BIN_DIR_TCL" OPTS=msvcrt,static,threads install
chmod -R a+rwx REPLACEME_SV_TOP_BIN_DIR_TCL
cd ../..
cd REPLACEME_SV_TK_DIR/win
nmake -f makefile.vc MACHINE=AMD64 TCLDIR="REPLACEME_SV_TOP_SRC_DIR_TCL" INSTALLDIR="REPLACEME_SV_TOP_BIN_DIR_TK" OPTS=msvcrt,static,threads hose
nmake -f makefile.vc MACHINE=AMD64 TCLDIR="REPLACEME_SV_TOP_SRC_DIR_TCL" INSTALLDIR="REPLACEME_SV_TOP_BIN_DIR_TK" OPTS=msvcrt,static,threads release
nmake -f makefile.vc MACHINE=AMD64 TCLDIR="REPLACEME_SV_TOP_SRC_DIR_TCL" INSTALLDIR="REPLACEME_SV_TOP_BIN_DIR_TK" OPTS=msvcrt,static,threads install
chmod -R a+rwx REPLACEME_SV_TOP_BIN_DIR_TK
nmake -f makefile.vc MACHINE=AMD64 TCLDIR="REPLACEME_SV_TOP_SRC_DIR_TCL" INSTALLDIR= OPTS=msvcrt,static,threads hose
cd ../..
cd REPLACEME_SV_TCL_DIR/win
nmake -f makefile.vc MACHINE=AMD64 INSTALLDIR="REPLACEME_SV_TOP_BIN_DIR_TCL" OPTS=msvcrt,static,threads hose

cd ../../BuildHelpers

export SV_TOPLEVEL_BINDIR_CYGWIN=`cygpath "REPLACEME_SV_TOP_BIN_DIR_TCL"`

cd ../tcllib-1.17
$SV_TOPLEVEL_BINDIR_CYGWIN/bin/REPLACEME_SV_TCLSH_EXECUTABLE installer.tcl
cd ../BuildHelpers
chmod -R a+rx $SV_TOPLEVEL_BINDIR_CYGWIN

cd ../tklib-0.6
$SV_TOPLEVEL_BINDIR_CYGWIN/bin/REPLACEME_SV_TCLSH_EXECUTABLE installer.tcl
cd ../BuildHelpers
chmod -R a+rx $SV_TOPLEVEL_BINDIR_CYGWIN

