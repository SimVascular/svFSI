## Compile `svFSI` on OSX

Follow the steps below to install `svFSI` on Mac.

1. Install Xcode command line tools using the command line below:

```bash
xcode-select install
```

2. Check if the default terminal is csh or bash. You may open the Terminal/Preferences and change the default shell to bash (/bin/bash). This is recommended as we may add some environmental variables later.

3. Install [Homebrew](https://docs.brew.sh/Installation)

4. Check if `/usr/local/bin` appears before `/usr/bin` in the `PATH` environmental variable.

```bash
echo $PATH
```

If not, open `~/.bash_profile` and add the following line:

```bash
export PATH="/usr/local/bin:$PATH"
```

5. Install the following packages using Homebrew

```bash
brew install gcc
brew install mpich
brew install cmake
```

Verify if MPI compilers are pointing to the correct compilers as,

```bash
which mpicc
which mpicxx
which mpif90
```

The CC and CXX compilers should point to the default `clang` compilers, while the FC/F90 compiler should point to `gfortran`. Set the environmental variable `LIBRARY_PATH` to the location of MPI_LIBRARY (usually `/usr/local/lib`) as,

```bash
export LIBRARY_PATH="/usr/local/lib:$LIBRARY_PATH"
```

in `~/.bash_profile`.

6. Mac OS usually comes with a default `blas` and `lapack` packages. If you prefer to use `OpenBLAS`, you may install it using Homebrew as,

```bash
brew install openblas
brew install lapack
```

You should also set the below environmental variables in `~/.bash_profile` as,

```bash
BLASDIR="/usr/local/Cellar/openblas/0.3.15_1"
export LIBRARY_PATH="$BLASDIR/lib:$LIBRARY_PATH"
```

Please make sure your to update the `BLASDIR` location depending on the version of OpenBLAS you install.

7. Download `svFSI` from GitHub:

```bash
mkdir svFSI-git
cd svFSI-git
git clone git@github.com:SimVascular/svFSI.git
mv svFSI/  src
```

8. Create a separate `build` folder and compile using the following commands:

```bash
mkdir build
cd build
ccmake ../
make
```

Enjoy!

