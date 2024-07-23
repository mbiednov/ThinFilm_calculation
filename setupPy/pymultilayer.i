%module pymultilayer

%{
    #define SWIG_FILE_WITH_INIT
    #include "pymultilayer.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%numpy_typemaps(bool, NPY_UINT, int)

%apply (double *IN_ARRAY1, int DIM1) {(double *n, int Nx1)}
%apply (double *IN_ARRAY1, int DIM1) {(double *k, int Nx2)}
%apply (double *IN_ARRAY1, int DIM1) {(double *thickness, int Nx3)}

%include "pymultilayer.h"