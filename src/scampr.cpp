#include <TMB.hpp>
#include "init.h" // for R CMD check: R_registerRoutines, R_useDynamicSymbols

template<class Type>
Type objective_function<Type>::operator() ()
{

  using namespace Eigen;

  // Utilities
  DATA_IVECTOR(bf_per_res);                 // Integer vector describing the number of basis functions within each spatial resolution
  DATA_INTEGER(mod_type);                   // Integer = 0 for no latent effects; = 1 for VA LGCP; = 2 for Laplace LGCP
  DATA_INTEGER(bf_type);                    // Integer = 0 for sparse; = 1 for dense basis function matrices
  DATA_INTEGER(data_type);                  // Integer = 0 for PO data; = 1 for PA data; = 2 for combined data

  // Parameters
  PARAMETER_VECTOR(fixed);                  // Fixed effects coefficients
  PARAMETER_VECTOR(random);                 // Random effect coefficients (for LAPLACE)
  PARAMETER_VECTOR(log_variance_component); // log(Variance) of each k random coefficient for VA; OR log(Standard Dev.) of random coefficient within spatial res. for LAPLACE - length l

  // Create enumeration of the data_type for switch function
  enum data_type { PO, PA, IDM };
  // Create enumeration of the mod_type for switch function
  enum mod_type { ipp, variational, laplace };
  // Create enumeration of the mod_type for switch function
  enum bf_type { sparse, dense };

  // Initialise the negative log-likelihood to be minimised (negative sum of the LL components)
  Type nll = 0.0;
  // Initialise the random effect contribution to the log-likelihood (as this is shared by each data type)
  Type LL_random = 0.0;
  // loop control variable for multi-resolution random effects
  int LCV = 0;

  switch (data_type){
  case PO:
  {
    #include "scampr_ppm.h"

    // Combine for the negative loglikelihood
    nll = -(LL_PO_pres + LL_PO_quad + LL_random);
    //REPORT(LL_PO_pres);
    //REPORT(LL_PO_quad);
    //REPORT(LL_random);

    break;
  }
  case PA:
  {
    #include "scampr_binom_cloglog.h"

    // Random coefficient likelihood component //
    // this is included outside 'scampr_binom_cloglog2.h' so it does not occur twice in the Integrated Data Model

    Type mu = 0.0;                                      // initialise the random coefficient mean
    vector<Type> PriorSD = exp(log_variance_component); // in the Laplace case this is exponentiated to get the std. dev. of random coeffs.

    for (int l = 0; l < bf_per_res.size(); ++l) {                   // loop through basis function resolutions
      for (int k = 0; k < bf_per_res(l); ++k) {                     // loop through random coefficients within resolution l
        LL_random += dnorm(random(k + LCV), mu, PriorSD(l), true);  // add the normal pdf contribution to the log-likelihood
      }
      LCV = LCV + bf_per_res(l);                                    // adjust the loop-control-variable so that next set of K_l coefficients are used
    }

    vector<Type> PriorVar = PriorSD * PriorSD;  // calculate the variance to match previous output
    // Add variance parameter to auto-diff reporting
    ADREPORT(PriorVar);

    // Combine for the loglikelihood
    nll = -(LL_PA + LL_random);
    //REPORT(LL_PA);
    //REPORT(LL_random);

    break;
  }
  case IDM:
  {
    #include "scampr_ppm.h"
    #include "scampr_binom_cloglog.h"

    // Combine for the loglikelihood
    nll = -(LL_PO_pres + LL_PO_quad + LL_PA + LL_random);
    //REPORT(LL_PO_pres);
    //REPORT(LL_PO_quad);
    //REPORT(LL_PA);
    //REPORT(LL_random);

    break;
  }
  default:
    error("Data type not recognised");
  } // end data type switch
  return nll;
}
