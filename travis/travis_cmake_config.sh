### install more recent version of CMake for Ubuntu 14.04
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
  wget http://simvascular.stanford.edu/downloads/public/open_source/linux/cmake/cmake-3.6.1-Linux-x86_64.sh
  chmod a+rx ./cmake-3.6.1-Linux-x86_64.sh
  sudo mkdir -p /usr/local/package/cmake-3.6.1
  sudo ./cmake-3.6.1-Linux-x86_64.sh --prefix=/usr/local/package/cmake-3.6.1 --skip-license
  sudo ln -s /usr/local/package/cmake-3.6.1/bin/ccmake    /usr/local/bin/ccmake
  sudo ln -s /usr/local/package/cmake-3.6.1/bin/cmake     /usr/local/bin/cmake
  sudo ln -s /usr/local/package/cmake-3.6.1/bin/cmake-gui /usr/local/bin/cmake-gui
  sudo ln -s /usr/local/package/cmake-3.6.1/bin/cpack     /usr/local/bin/cpack
  sudo ln -s /usr/local/package/cmake-3.6.1/bin/ctest     /usr/local/bin/ctest
fi

#compilers
if [[ "$TRAVIS_OS_NAME" == "linux" ]]
then
  export CC="gcc"
  export CXX="g++"
elif [[ "$TRAVIS_OS_NAME" == "osx" ]]
then
  export CC="clang"
  export CXX="clang++"
fi

pushd $BUILD_DIR

echo SV_EXTERNALS_ARGS: $SV_EXTERNALS_ARGS
#cmake
export REPLACEME_SV_CMAKE_CMD="cmake"
export REPLACEME_SV_CMAKE_GENERATOR="Unix Makefiles"
export REPLACEME_SV_CMAKE_BUILD_TYPE="RelWithDebInfo"
export REPLACEME_SV_MAKE_CMD="make -j8"
export REPLACEME_SV_TOP_SRC_DIR_SV=$SV_TOP_DIR

"$REPLACEME_SV_CMAKE_CMD" \
\
   -G "$REPLACEME_SV_CMAKE_GENERATOR" \
\
   -DCMAKE_BUILD_TYPE="$REPLACEME_SV_CMAKE_BUILD_TYPE" \
   -DBUILD_SHARED_LIBS=OFF \
   -DBUILD_TESTING=OFF \
   -DSV_EXTERNALS_ADDITIONAL_CMAKE_ARGS="$SV_EXTERNALS_ARGS" \
\
 "$REPLACEME_SV_TOP_SRC_DIR_SV"

popd
