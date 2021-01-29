#include <TMB.hpp>
#include "init.h" // for R CMD check: R_registerRoutines, R_useDynamicSymbols

template<class Type>
  Type objective_function<Type>::operator() ()
{

    using namespace Eigen;

    DATA_STRING(bf_matrix_type);
    if (bf_matrix_type == "dense") {
      #include "popa_dense.h"
    } else if (bf_matrix_type == "sparse") {
      #include "popa_sparse.h"
    } else {
      error ("bf_matrix_type must be 'dense' or 'sparse'");
    }
}
