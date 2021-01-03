  // Presence-only Data
  DATA_MATRIX(X_PO_pres);                   // Fixed effects n_{PO} x p matrix - at the presence pts
  DATA_MATRIX(B_PO_pres);                   // Biasing effects n_{PO} x b matrix - at the presence pts
  DATA_MATRIX(Z_PO_pres);            // Basis functions n_{PO} x k matrix - at the presence pts
  DATA_MATRIX(X_PO_quad);                   // Fixed effects m x p matrix - at the quadrature pts
  DATA_MATRIX(B_PO_quad);                   // Biasing effects m x b matrix - at the quadrature pts
  DATA_MATRIX(Z_PO_quad);            // Basis functions m x k matrix - at the quadrature pts
  DATA_VECTOR(quad_size);                   // Vector describing the quadrature sizes - length m
  // Presence-Absence Data
  DATA_MATRIX(X_PA);                        // Fixed effects n_{PA} x p matrix - at the sites
  DATA_MATRIX(Z_PA);                 // Basis functions n_{PA} x k matrix - at the sites
  DATA_VECTOR(Y);                           // Binary vector describing presence/absence - length n_{PA}
  // Utilities
  DATA_IVECTOR(bf_per_res);                 // Integer vector describing the number of basis functions within each spatial resolution
  DATA_INTEGER(mod_type);                   // Integer = 0 for no latent effects; = 1 for VA LGCP; = 2 for Laplace LGCP
  DATA_INTEGER(data_type);                  // Integer = 0 for PO data; = 1 for PA data; = 2 for combined data
  // Parameters
  PARAMETER_VECTOR(fixed);                  // Fixed effects coefficients
  PARAMETER_VECTOR(bias);                   // Biasing effect coefficients
  // COULD WE COMBINE THESE? PROBLEM AS DIMENSION REQ. FOR X * fixed WOULD CHANGE FOR PO & POPA
  PARAMETER_VECTOR(random);                 // Random effect coefficients (for LAPLACE) OR coefficient means (for VARIATIONAL)
  PARAMETER_VECTOR(log_variance_component); // log(Standard Dev.) of Random coefficient within spatial res. (for LAPLACE - length l) OR log(Variances) of VA Random effects (for VARIATIONAL  - length k)

  // Create enumeration of the data_type for switch function
  enum data_type { PO, PA, POPA };
  
  // Create enumeration of the mod_type for switch function
  enum mod_type { ipp, variational, laplace };
  
  // Initialise the various log-likelihood components
  Type LL_PO_pres = 0.0;
  Type LL_PO_quad = 0.0;
  Type LL_PA = 0.0;
  Type LL_random = 0.0;
  
  // Initialise the negative log-likelihood to be minimised 
  Type nll = 0.0;
  
  // loop control variable for multi-resolution models
  int L = 0;
  
  // Switch for different data types being used
  switch(data_type){
  case PO:
  {
    // Shared by all model types:
    vector<Type> Xfixed_PO_pres = X_PO_pres * fixed;
    vector<Type> Bbias_PO_pres = B_PO_pres * bias;
    vector<Type> Xfixed_PO_quad = X_PO_quad * fixed;
    vector<Type> Bbias_PO_quad = B_PO_quad * bias;
    vector<Type> PRES = Xfixed_PO_pres + Bbias_PO_pres;
    vector<Type> QUAD = Xfixed_PO_quad + Bbias_PO_quad;
    
    switch (mod_type){
    case ipp:
    {
      LL_PO_pres += PRES.sum();
      LL_PO_quad -= (quad_size * exp(QUAD)).sum();
      break;
    }
    case variational:
    {
      // log-likelihood components at the presence pts and quadrature //
      vector<Type> PosteriorVar = exp(log_variance_component);
      vector<Type> Zrandom_mean_pres = Z_PO_pres * random;
      vector<Type> Zrandom_mean_quad = Z_PO_quad * random;
      matrix<Type> Z2 = Z_PO_quad.cwiseProduct(Z_PO_quad);
      vector<Type> ZsigZ = Z2 * PosteriorVar;
      LL_PO_pres += PRES.sum() + Zrandom_mean_pres.sum();
      LL_PO_quad -= (quad_size * exp(QUAD + Zrandom_mean_quad + (0.5 * ZsigZ))).sum();
      
      // Random coefficient likelihood component //
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
      LL_random = -0.5 * DKL.sum();
      ADREPORT(PosteriorVar);
      ADREPORT(PriorVar);
      break;
    }
    case laplace:
    {
      // log-likelihood components at the presence pts and quadrature //
      
      vector<Type> Zrandom_pres = Z_PO_pres * random;
      vector<Type> Zrandom_quad = Z_PO_quad * random;
      LL_PO_pres += PRES.sum() + Zrandom_pres.sum();
      LL_PO_quad -= (quad_size * exp(QUAD + Zrandom_quad)).sum();
      
      // Random coefficient likelihood component //
      
      Type mu = 0.0; // initialise the random coefficient mean
      vector<Type> PriorSD = exp(log_variance_component);
      L = 0;
      for (int l = 0; l < bf_per_res.size(); ++l) {
        for (int k = 0; k < bf_per_res(l); ++k) {
          LL_random += dnorm(random(k + L), mu, PriorSD(l), true);
        }
        L = L + bf_per_res(l);
      }
      vector<Type> PriorVar = PriorSD * PriorSD;
      ADREPORT(PriorVar);
      break;
    }
    default:
      error("Model Type not recognised");
    } // end model type switch
    break;
  }
  case PA:
  {
    // Matrix operations shared by all model types
    vector<Type> Xfixed_PA = X_PA * fixed;
    vector<Type> prob_abs(Xfixed_PA.size());
      
    switch (mod_type){
    case ipp:
    {
      // probability of absence at the sites  //
      prob_abs = exp(-exp(Xfixed_PA));
      break;
    }
    case variational: // NO CLOSED FORM SOLUTION! SAME AS LAPLACE FOR NOW
    {
      // probability of absence at the sites  //
      vector<Type> Zrandom_PA = Z_PA * random;
      prob_abs = exp(-exp(Xfixed_PA + Zrandom_PA));
      
      // Random coefficient likelihood component //
      
      Type mu = 0.0; // initialise the random coefficient mean
      vector<Type> PriorSD = exp(log_variance_component);
      L = 0;
      for (int l = 0; l < bf_per_res.size(); ++l) {
        for (int k = 0; k < bf_per_res(l); ++k) {
          LL_random += dnorm(random(k + L), mu, PriorSD(l), true);
        }
        L = L + bf_per_res(l);
      }
      
      vector<Type> PriorVar = PriorSD * PriorSD;
      ADREPORT(PriorVar);
      break;
    }
    case laplace:
    {
      // probability of absence at the sites  //
      vector<Type> Zrandom_PA = Z_PA * random;
      prob_abs = exp(-exp(Xfixed_PA + Zrandom_PA));
      
      // Random coefficient likelihood component //
      
      Type mu = 0.0; // initialise the random coefficient mean
      vector<Type> PriorSD = exp(log_variance_component);
      L = 0;
      for (int l = 0; l < bf_per_res.size(); ++l) {
        for (int k = 0; k < bf_per_res(l); ++k) {
          LL_random += dnorm(random(k + L), mu, PriorSD(l), true);
        }
        L = L + bf_per_res(l);
      }
      
      vector<Type> PriorVar = PriorSD * PriorSD;
      ADREPORT(PriorVar);
      break;
    }
    default:
      error("Model Type not recognised");
    } // end model type switch
    
    // log-likelihood component at the presence/absence sites //
    
    Type Size = 1;
    for (int i = 0; i < X_PA.rows(); i++) {
      LL_PA += dbinom(Y(i), Size, 1 - prob_abs(i), true);
    }
    break;
  }
  case POPA:
  {
    // Shared by all model types:
    vector<Type> Xfixed_PO_pres = X_PO_pres * fixed;
    vector<Type> Bbias_PO_pres = B_PO_pres * bias;
    vector<Type> Xfixed_PO_quad = X_PO_quad * fixed;
    vector<Type> Bbias_PO_quad = B_PO_quad * bias;
    vector<Type> PRES = Xfixed_PO_pres + Bbias_PO_pres;
    vector<Type> QUAD = Xfixed_PO_quad + Bbias_PO_quad;
    vector<Type> Xfixed_PA = X_PA * fixed;
    vector<Type> prob_abs(Xfixed_PA.size());
    
    switch (mod_type){
    case ipp:
    {
      // log-likelihood component at the presence and quad points //
      LL_PO_pres += PRES.sum();
      LL_PO_quad -= (quad_size * exp(QUAD)).sum();

      // probability of absence at the sites  //
      prob_abs = exp(-exp(Xfixed_PA));
      break;
    }
    case variational:  // NO CLOSED FORM SOLUTION! SAME AS LAPLACE FOR NOW
    {
      // log-likelihood components at the presence pts and quadrature //
      
      vector<Type> Zrandom_pres = Z_PO_pres * random;
      vector<Type> Zrandom_quad = Z_PO_quad * random;
      LL_PO_pres += PRES.sum() + Zrandom_pres.sum();
      LL_PO_quad -= (quad_size * exp(QUAD + Zrandom_quad)).sum();
      
      // probability of absence at the sites  //
      
      vector<Type> Zrandom_PA = Z_PA * random;
      prob_abs = exp(-exp(Xfixed_PA + Zrandom_PA));
      
      // Random coefficient likelihood component //
      
      Type mu = 0.0; // initialise the random coefficient mean
      vector<Type> PriorSD = exp(log_variance_component);
      L = 0;
      for (int l = 0; l < bf_per_res.size(); ++l) {
        for (int k = 0; k < bf_per_res(l); ++k) {
          LL_random += dnorm(random(k + L), mu, PriorSD(l), true);
        }
        L = L + bf_per_res(l);
      }
      vector<Type> PriorVar = PriorSD * PriorSD;
      ADREPORT(PriorVar);
      break;
    }
    case laplace:
    {
      // log-likelihood components at the presence pts and quadrature //
      
      vector<Type> Zrandom_pres = Z_PO_pres * random;
      vector<Type> Zrandom_quad = Z_PO_quad * random;
      LL_PO_pres += PRES.sum() + Zrandom_pres.sum();
      LL_PO_quad -= (quad_size * exp(QUAD + Zrandom_quad)).sum();
      
      // probability of absence at the sites  //
      
      vector<Type> Zrandom_PA = Z_PA * random;
      prob_abs = exp(-exp(Xfixed_PA + Zrandom_PA));
      
      // Random coefficient likelihood component //
      
      Type mu = 0.0; // initialise the random coefficient mean
      vector<Type> PriorSD = exp(log_variance_component);
      L = 0;
      for (int l = 0; l < bf_per_res.size(); ++l) {
        for (int k = 0; k < bf_per_res(l); ++k) {
          LL_random += dnorm(random(k + L), mu, PriorSD(l), true);
        }
        L = L + bf_per_res(l);
      }
      vector<Type> PriorVar = PriorSD * PriorSD;
      ADREPORT(PriorVar);
      break;
    }
    default:
      error("Model Type not recognised");
    } // end model type switch
    
    // log-likelihood component at the presence/absence sites //
    
    Type Size = 1;
    for (int i = 0; i < X_PA.rows(); i++) {
      LL_PA += dbinom(Y(i), Size, 1 - prob_abs(i), true);
    }
    break;
  }
  default:
    error("Data Type not recognised");
  } // end data type switch
  
  // Combine for the loglikelihood
  nll = -(LL_PO_pres + LL_PO_quad + LL_PA + LL_random);
  REPORT(LL_PO_pres);
  REPORT(LL_PO_quad);
  REPORT(LL_PA);
  REPORT(LL_random);
  return nll;