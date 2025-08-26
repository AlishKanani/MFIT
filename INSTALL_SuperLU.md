# Installing SuperLU as Dynamic Library

<!-- note these instructions are valid for linux system -->

Note: below instructions are tested on linux system.

1. Download the SuperLU source code from [here](https://github.com/xiaoyeli/superlu).

```bash
wget https://github.com/xiaoyeli/superlu/archive/refs/tags/v7.0.0.zip
unzip v7.0.0.zip
cd superlu-7.0.0
```
Make sure you have `cmake` installed on your system. 

2. Create a build directory and run cmake.

```bash
mkdir build; cd build
cmake -DTPL_BLAS_LIBRARIES=/lib/x86_64-linux-gnu/libblas.so -DBUILD_SHARED_LIBS=TRUE ..
```

DBUILD_SHARED_LIBS=TRUE is used to build the SuperLU library as a dynamic library.

Replace `/lib/x86_64-linux-gnu/libblas.so` with the path to your dynamic library of BLAS. You can find the path using `ldconfig -p | grep blas`.

3. Build the SuperLU library.

```bash
make
make install
```

4. Add the path to the SuperLU library and blas library to your `LD_LIBRARY_PATH`.

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lib/x86_64-linux-gnu/
```
Replace `/usr/local/lib/` with the path to the directory where the SuperLU library is installed.
