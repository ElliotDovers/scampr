// Laplace Approx. of log-likelihood for presence-only data Cox Process
// Single or Multiple Resolution, non-spatially correlated random coefficients
// Sparse Basis Functions
#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() ()
{

  using namespace Eigen;

  DATA_MATRIX(X_pres);               // Fixed effect (N+M)xP matrix - at the pres & quad pts
  DATA_SPARSE_MATRIX(Z_pres);        // Basis functions (N+M)xK matrix - at the pres & quad pts
  DATA_MATRIX(X_quad);               // Fixed effect (N+M)xP matrix - at the pres & quad pts
  DATA_SPARSE_MATRIX(Z_quad);        // Basis functions (N+M)xK matrix - at the pres & quad pts
  DATA_VECTOR(quad_size);            // Vector describing the quadrature sizes - length N + M
  PARAMETER_VECTOR(fixed);           // Fixed effect coefficients
  PARAMETER_VECTOR(random);          // Random effect coefficients (for LAPLACE) | VARIATIONAL coefficient means
  PARAMETER_VECTOR(log_variance_component);      // LAPLACE: log(Standard Dev. of Random coefficient within spatial res. - length L | VARIATIONAL: Random effect VA variances - length K
  DATA_IVECTOR(bf_per_res);          // A utility integer vector describing the number of basis functions within each spatial res.
  DATA_INTEGER(mod_type);            // 0 for ipp, 1 for VA LGCP, 2 for Laplace LGCP

  Type ll1 = 0.0;   // initialise the log-likelihood component 1
  Type ll2 = 0.0;   // initialise the log-likelihood component 2
  Type ll3 = 0.0;   // initialise the log-likelihood component 3
  int L = 0;        // loop control variable for LGCP models

  // Create enumeration of the mod_type factor
  enum mod_type { ipp, variational, laplace };

  // Shared by all models:
  vector<Type> Xfixed_pres = X_pres * fixed;
  vector<Type> Xfixed_quad = X_quad * fixed;

  switch (mod_type){
  case ipp:
  {
    ll1 += Xfixed_pres.sum();
    ll2 -= (quad_size * exp(Xfixed_quad)).sum();
    break;
  }
  case variational:
  {
    // int K = random.size();
    vector<Type> PosteriorVar = exp(log_variance_component);
    vector<Type> Zrandom_mean_pres = Z_pres * random;
    vector<Type> Zrandom_mean_quad = Z_quad * random;
    SparseMatrix<Type> Z2= Z_quad.cwiseProduct(Z_quad);
    vector<Type> ZsigZ = Z2 * PosteriorVar;
    ll1 += Xfixed_pres.sum() + Zrandom_mean_pres.sum();
    ll2 -= (quad_size * exp(Xfixed_quad + Zrandom_mean_quad + (0.5 * ZsigZ))).sum();
    L = 0;
    vector<Type> PriorVar(bf_per_res.size());
    vector<Type> DKL(random.size());
    for (int l = 0; l < bf_per_res.size(); ++l) {
      vector<Type> tmp(bf_per_res(l));
      for (int k = 0; k < bf_per_res(l); ++k) {
        tmp(k) = PosteriorVar(k + L) + (random(k + L) * random(k + L));
      }
      PriorVar(l) = tmp.sum() / tmp.size();
      for (int k = 0; k < bf_per_res(l); ++k) {
        DKL(k + L) = log(PriorVar(l)) - log_variance_component(k + L);
      }
      L = L + bf_per_res(l);
    }
    ll3 = -0.5 * DKL.sum();
    ADREPORT(PosteriorVar);
    ADREPORT(PriorVar);
    break;
  }
  case laplace:
  {
    Type mu = 0.0;    // initialise the random coefficient mean
    vector<Type> PriorSD = exp(log_variance_component);
    vector<Type> Zrandom_pres = Z_pres * random;
    vector<Type> Zrandom_quad = Z_quad * random;
    ll1 += Xfixed_pres.sum() + Zrandom_pres.sum();
    ll2 -= (quad_size * exp(Xfixed_quad + Zrandom_quad)).sum();
    L = 0;
    for (int l = 0; l < bf_per_res.size(); ++l) {
      for (int k = 0; k < bf_per_res(l); ++k) {
        ll3 += dnorm(random(k + L), mu, PriorSD(l), true);
      }
      L = L + bf_per_res(l);
    }
    vector<Type> PriorVar = PriorSD * PriorSD;
    ADREPORT(PriorVar);
    break;
  }
  default:
    error("Model Type not recognised");
  } // end switch

  // Combine for the loglikelihood
  Type nll = -(ll1 + ll2 + ll3);
  REPORT(ll1);
  REPORT(ll2);
  REPORT(ll3);

  return nll;
  }
