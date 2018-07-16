# QDLDL
A free LDL factorisation routine for quasi-definite linear systems: `Ax=b`

[![Build Status](https://travis-ci.org/oxfordcontrol/qdldl.svg?branch=master)](https://travis-ci.org/oxfordcontrol/qdldl)
[![Build status](https://ci.appveyor.com/api/projects/status/ns4br7v6y3i5stai/branch/master?svg=true)](https://ci.appveyor.com/project/bstellato/qdldl-8q1mv/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/oxfordcontrol/qdldl/badge.svg)](https://coveralls.io/github/oxfordcontrol/qdldl)

QDLDL is a new free linear system solver for quasi-definite linear systems.


## Getting started
To start using QDLDL just clone the repository

```bash
git clone https://github.com/oxfordcontrol/qdldl.git
```

### Build

To build it, you need to install [cmake](https://cmake.org/) and run

```bash
mkdir build/
cd build/
cmake ..
make
```

This will generate an `out/` folder where you will have:

- `qdldl_example`: a **code example** from [`examples/c/example.c`](./examples/c/example.c)
- `libqdldlstatic`: a static library
- `libqdldl`: a dynamic library

**N.B.** All files will have the extension relative to the operating system used.


### Install/Uninstall

To install (uninstall) the libraries and headers you can simply run `make install` (`make uninstall`) after running the `cmake` command above.


## Calling QDLDL

### Including the header
QDLDL uses its internal types for integer floats and booleans. If you want to overwrite them with your worn you need to redefine them before including the library
```c
#typedef mybool QDLDL_bool;
#typedef myint QDLDL_int;
#typedef myfloat QDLDL_float;

/* Need to specify this to avoid types redefinition */
#define QDLDL_TYPES_DEFINED

#include "qdldl.h"
```

### Main API

The QDLDL API consists in 5 functions documented in `include/qdldl.h`.
For more details see the example in [`examples/c/example.c`](./examples/c/example.c).

**N.B.** There is **no memory allocation** performed in these routines. The user is assumed to have the working vectors already allocated.

Here is a brief summary.

* `QDLDL_etree`: compute the elimination tree for the factorization `A = LDL'`
* `QDLDL_factor`: return the factors `L`, `D` and `Dinv = 1./D`
* `QDLDL_solve`: solve the linear system `LDL'x = b`
* `QDLDL_Lsolve`: solve `(L + I)x = b`
* `QDLDL_Ltsolve`: solve `(L + I)'x = b`


## Linking QDLDL

### Basic Example
A basic example appears in [`examples/c/example.c`](./examples/c/example.c) and is compiled using cmake and the `CMakeLists.txt` file in the root folder


### Including in a cmake project

You can include QDLDL in a cmake project `foo` by adding the subdirectory as
```
# Add project
add_subdirectory(qdldl EXCLUDE_FROM_ALL)
```

QDLDL can be linked using a static or dynamic linker
```
# Link static library
target_link_libraries (foo qdldlstatic)

# Link shared library
target_link_libraries (foo qdldl)
```
for dynamic linking the shared library should be available in your path.


## The algorithm

The algorithm is an independent implementation of the elimination tree and factoring procedures proposed in

> T. A Davis. [Algorithm 849: a concise sparse cholesky factorization package](https://dl.acm.org/citation.cfm?id=1114277). ACM Trans. Math. Softw., 31(4):587–591, 2005.


## Credits

- [Paul Goulart](http://users.ox.ac.uk/~engs1373/): main development
- [Bartolomeo Stellato](https://stellato.io/): code refactoring and testing
- [Goran Banjac](http://control.ee.ethz.ch/~gbanjac/)







