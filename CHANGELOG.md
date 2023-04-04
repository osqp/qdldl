Version 0.1.7 (4 April 2023)
---------------------------------
*   Using the correct CMake variable for longs during compilation
*   Using the correct CMake variable for longs during CI Testing 

Version 0.1.6 (20 June 2022)
---------------------------------
*   Export the version number in CMake and in preprocessor macros.
*   Add CMake options to enable/disable building the shared/static libraries and demo executable.
*   Add symbol visibility information to the public API and build a Windows import library.
*   Rename all CMake options to have a `QDLDL_` prefix, e.g.
    * `DFLOAT` -> `QDLDL_FLOAT`
    * `DLONG` -> `QDLDL_LONG`
    * `UNITTESTS` -> `QDLDL_UNITTESTS`


Version 0.1.5 (7 August 2020)
---------------------------------
*   Reduce the amount of memory access in each iteration of the solver loops.


Version 0.1.4 (6 September 2019)
---------------------------------
*   Fix to cmake when building `qdldl::qdldlstatic`.
*   Fix to overflow issue when factoring very large matrices.


Version 0.1.3 (11 September 2018)
----------------------------------
*   Julia implementation supports logical factorisation.
*   Changed `QDLDL_bool` to be unsigned char everywhere, except for the Julia examples where it is now treated as `Uint8`.


Version 0.1.2 (23 July 2018)
-----------------------------
*   Various cmake improvements.
*   Added pure Julia implementation.


Version 0.1.1 (19 July 2018)
-----------------------------
*   Fixed behaviour when data in A does not appear
    sequentially within each column.
*   Additional unit tests for non-sequential columns.
*   Types can be defined through cmake.


Version 0.1.0 (16 July 2018)
-----------------------------
*   Initial release
