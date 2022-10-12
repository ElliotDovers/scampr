  // Point Process Model //
  // Either an inhomogeneous Poisson Process or log-Gaussian Cox Process (the latter using a Laplace or variational approx.)
  
  // Presence-only Data
  DATA_MATRIX(X_PO_pres);                   // Fixed effects n_{PO} x p matrix - at the presence pts
  DATA_MATRIX(X_PO_quad);                   // Fixed effects q x p matrix - at the quadrature pts
  DATA_VECTOR(quad_size);                   // Vector describing the quadrature sizes - length m
  
  // Utility specific to PPM
  DATA_INTEGER(bias_type);                   // Integer = 0 for none (or included in fixed effects); = 1 for biasing covariates; = 2 for secondary latent field (IDM only)
  // Create enumeration of the bias_type for switch function
  enum bias_type { none, covars, latent, new_latent};

  // Initialise the various log-likelihood components (these are simply added to or subtracted from later)
  Type LL_PO_pres = 0.0;
  Type LL_PO_quad = 0.0;
  
  // Shared by all model types:
  vector<Type> Xfixed_PO_pres = X_PO_pres * fixed;  // fixed effects at presence pts. (including intercept depending on formula)
  vector<Type> Xfixed_PO_quad = X_PO_quad * fixed;  // fixed effects at quadrature pts. (including intercept depending on formula)
  
  // Initialise the vectors resulting from the switch (will add in bias as needed)
  vector<Type> eta_fixed_pres = Xfixed_PO_pres; // at presence pts.
  vector<Type> eta_fixed_quad = Xfixed_PO_quad; // at quadrature pts.
  
  // Calculate the bias components when included as covariates (separate to fixed effects)
  if (bias_type == covars) {
    // Data specific to PPM with biasing covariates
    DATA_MATRIX(B_PO_pres);                           // Biasing effects n_{PO} x b matrix - at the presence pts
    DATA_MATRIX(B_PO_quad);                   // Biasing effects q x b matrix - at the quadrature pts
    // Parameters specific to PPM with biasing covariates
    PARAMETER_VECTOR(bias);                   // Biasing effect coefficients
    
    vector<Type> Bbias_PO_pres = B_PO_pres * bias;    // biasing effects at presence pts. (including intercept depending on formula)
    vector<Type> Bbias_PO_quad = B_PO_quad * bias;    // biasing effects at quadrature pts. (including intercept depending on formula)
    
    // add bias component to linear predictor at presence and quadrature pts.
    eta_fixed_pres = eta_fixed_pres + Bbias_PO_pres;
    eta_fixed_quad = eta_fixed_quad + Bbias_PO_quad;
  }

  switch (mod_type){
  case ipp:
  {
    LL_PO_pres += eta_fixed_pres.sum();                     // sum of log-intensity at the presence pts.
    LL_PO_quad -= (quad_size * exp(eta_fixed_quad)).sum();  // approx. of spatial integral of the intensity
    break;
  }
  case variational:
  {
    // log-likelihood components at the presence and quadrature points //
    
    vector<Type> PosteriorVar = exp(log_variance_component); // in the VA case this is exponentiated to get the variance on the random coeffs.
    
    // Initialise the vectors resulting from the switch
    vector<Type> Zrandom_mean_pres; // latent field approx. at the presence points
    vector<Type> Zrandom_mean_quad; // latent field approx. at the quadrature points
    vector<Type> ZsigZ;             // Z Sigma t(Z) at the quadrature points (note this is quick due to diagonal Sigma)
    
    // Calculate the random components based on switch
    switch (bf_type){
    case sparse:
    {
      DATA_SPARSE_MATRIX(Z_PO_pres);            // Basis functions n_{PO} x k matrix - at the presence pts
      DATA_SPARSE_MATRIX(Z_PO_quad);            // Basis functions m x k matrix - at the quadrature pts
      
      Zrandom_mean_pres = Z_PO_pres * random;                     // random effects at presence pts.
      Zrandom_mean_quad = Z_PO_quad * random;                     // fixed effects at quadrature pts.
      SparseMatrix<Type> Z2 = Z_PO_quad.cwiseProduct(Z_PO_quad);  // using Eigen to quickly calc. element-wise squaring
      ZsigZ = Z2 * PosteriorVar;                                  // correction in MGF for the intensity
      break;
    }
    case dense:
    {
      DATA_MATRIX(Z_PO_pres);            // Basis functions n_{PO} x k matrix - at the presence pts
      DATA_MATRIX(Z_PO_quad);            // Basis functions m x k matrix - at the quadrature pts
      
      Zrandom_mean_pres = Z_PO_pres * random;               // random effects at presence pts.
      Zrandom_mean_quad = Z_PO_quad * random;               // fixed effects at quadrature pts.
      matrix<Type> Z2 = Z_PO_quad.cwiseProduct(Z_PO_quad);  // using Eigen to quickly calc. element-wise squaring
      ZsigZ = Z2 * PosteriorVar;                            // correction in MGF for the intensity
      break;
    }  default:
      error("Basis function matrix type not recognised");
    } // end bf_type switch
    

    LL_PO_pres += eta_fixed_pres.sum() + Zrandom_mean_pres.sum();                               // sum of log-intensity at the presence pts.
    LL_PO_quad -= (quad_size * exp(eta_fixed_quad + Zrandom_mean_quad + (0.5 * ZsigZ))).sum();  // approx. of spatial integral of the MGF of intensity

    // Random coefficient likelihood component //
    
    vector<Type> PriorVar(bf_per_res.size());                                 // initialise the prior variance on the random coefficients
    vector<Type> DKL(random.size());                                          // initialise the Kullback-Leibler div. between va density and prior density of random coeffs.
    for (int l = 0; l < bf_per_res.size(); ++l) {                             // loop through basis function resolutions
      vector<Type> tmp(bf_per_res(l));                                        // initialise a temporary vector to calculate the MLE of within level prior variance
      for (int k = 0; k < bf_per_res(l); ++k) {                               // loop through each random coefficient k within resolution, l
        tmp(k) = PosteriorVar(k + LCV) + (random(k + LCV) * random(k + LCV)); // calc. sigma_k^2 + mu_k^2
      }
      PriorVar(l) = tmp.sum() / tmp.size();                                   // calc. MLE of prior variance within resolution level, l
      for (int k = 0; k < bf_per_res(l); ++k) {                               // loop through each random coefficient k within resolution, l
        DKL(k + LCV) = log(PriorVar(l)) - log_variance_component(k + LCV);    // calc. Kullback-Leibler using MLE of prior variance
      }
      LCV = LCV + bf_per_res(l);                                              // adjust the loop-control-variable so that next set of K_l coefficients are used
    }
    
    LL_random = -0.5 * DKL.sum(); // random component of log-likelihood
    
    // Add variance parameters to auto-diff reporting
    ADREPORT(PosteriorVar);
    ADREPORT(PriorVar);
    break;
  }
  case laplace:
  {
    // log-likelihood components at the presence pts and quadrature //
    
    // Initialise the vectors resulting from the switch
    vector<Type> Zrandom_pres; // latent field approx. at the presence points
    vector<Type> Zrandom_quad; // latent field approx. at the quadrature points
    
    // Calculate the random components based on switch
    switch (bf_type){
    case sparse:
    {
      DATA_SPARSE_MATRIX(Z_PO_pres);      // Basis functions n_{PO} x k matrix - at the presence pts
      DATA_SPARSE_MATRIX(Z_PO_quad);      // Basis functions m x k matrix - at the quadrature pts
    
      Zrandom_pres = Z_PO_pres * random;  // random effects at presence pts.
      Zrandom_quad = Z_PO_quad * random;  // random effects at quadrature pts.
      
      // Calculate the random biasing components based on switch
      switch (bias_type){
      case latent:
      {
        // Parameters specific to IDM accounting for bias with a secondary latent field
        PARAMETER_VECTOR(bias); // Biasing coefficients (are random in this instance)
        
        vector<Type> Zrandom_bias_pres = Z_PO_pres * bias; // random biasing field at presence pts.
        vector<Type> Zrandom_bias_quad = Z_PO_quad * bias; // random biasing field at quadrature pts.
        
        // add on the random biasing fields
        Zrandom_pres = Zrandom_pres + Zrandom_bias_pres;
        Zrandom_quad = Zrandom_quad + Zrandom_bias_quad;
        break;
      }
      case new_latent:
      {
        // Parameters specific to IDM accounting for bias with a secondary latent field
        PARAMETER_VECTOR(bias); // Biasing coefficients (are random in this instance)
        
        DATA_SPARSE_MATRIX(Z2_PO_pres);        // Secondary basis functions n_{PO} x k_2 matrix - at the presence pts
        DATA_SPARSE_MATRIX(Z2_PO_quad);        // Secondary basis functions m x k_2 matrix - at the quadrature pts
        
        vector<Type> Zrandom_bias_pres = Z2_PO_pres * bias; // random biasing field at presence pts.
        vector<Type> Zrandom_bias_quad = Z2_PO_quad * bias; // random biasing field at quadrature pts.
        
        // add on the random biasing fields
        Zrandom_pres = Zrandom_pres + Zrandom_bias_pres;
        Zrandom_quad = Zrandom_quad + Zrandom_bias_quad;
        break;
      }
      } // end of basis_type switch
      
      break;
    }
    case dense:
    {
      DATA_MATRIX(Z_PO_pres);            // Basis functions n_{PO} x k matrix - at the presence pts
      DATA_MATRIX(Z_PO_quad);            // Basis functions m x k matrix - at the quadrature pts
    
      Zrandom_pres = Z_PO_pres * random;  // random effects at presence pts.
      Zrandom_quad = Z_PO_quad * random;  // random effects at quadrature pts.
      
      // Calculate the random biasing components based on switch
      switch (bias_type){
      case latent:
      {
        // Parameters specific to IDM accounting for bias with a secondary latent field
        PARAMETER_VECTOR(bias); // Biasing coefficients (are random in this instance)
        
        
        vector<Type> Zrandom_bias_pres = Z_PO_pres * bias; // random biasing field at presence pts.
        vector<Type> Zrandom_bias_quad = Z_PO_quad * bias; // random biasing field at quadrature pts.
        
        // add on the random biasing fields
        Zrandom_pres = Zrandom_pres + Zrandom_bias_pres;
        Zrandom_quad = Zrandom_quad + Zrandom_bias_quad;
        break;
      }
      case new_latent:
      {
        // Parameters specific to IDM accounting for bias with a secondary latent field
        PARAMETER_VECTOR(bias); // Biasing coefficients (are random in this instance)
        
        DATA_MATRIX(Z2_PO_pres);        // Secondary basis functions n_{PO} x k_2 matrix - at the presence pts
        DATA_MATRIX(Z2_PO_quad);        // Secondary basis functions m x k_2 matrix - at the quadrature pts
        
        vector<Type> Zrandom_bias_pres = Z2_PO_pres * bias; // random biasing field at presence pts.
        vector<Type> Zrandom_bias_quad = Z2_PO_quad * bias; // random biasing field at quadrature pts.
        
        // add on the random biasing fields
        Zrandom_pres = Zrandom_pres + Zrandom_bias_pres;
        Zrandom_quad = Zrandom_quad + Zrandom_bias_quad;
        break;
      }
      } // end of basis_type switch
    
      break;
    }
    default:
      error("Basis function matrix type not recognised");
    } // end bf_type switch

    LL_PO_pres += eta_fixed_pres.sum() + Zrandom_pres.sum();              // sum of log-intensity at the presence pts.
    LL_PO_quad -= (quad_size * exp(eta_fixed_quad + Zrandom_quad)).sum(); // approx. of spatial integral of the intensity

    // Random coefficient likelihood component //

    Type mu = 0.0;                                      // initialise the random coefficient mean
    vector<Type> PriorSD = exp(log_variance_component); // in the Laplace case this is exponentiated to get the std. dev. of random coefficients

    for (int l = 0; l < bf_per_res.size(); ++l) {                   // loop through basis function resolutions
      for (int k = 0; k < bf_per_res(l); ++k) {                     // loop through random coefficients within resolution l
        LL_random += dnorm(random(k + LCV), mu, PriorSD(l), true);  // add the normal pdf contribution to the log-likelihood
      }
      LCV = LCV + bf_per_res(l);                                    // adjust the loop-control-variable so that next set of K_l coefficients are used
    }
    
    vector<Type> PriorVar = PriorSD * PriorSD;  // calculate the variance to match that of the VA case
    // Add variance parameter to auto-diff reporting
    ADREPORT(PriorVar);
    
    // Add on additional random coefficient contributions to LL_random when using an IDM with secondary latent field
    switch (bias_type){
    case latent:
    {
      // Parameters specific to IDM accounting for bias with a secondary latent field
      PARAMETER_VECTOR(log_variance_component_bias); // log(Standard Dev.) of random coefficients within spatial res., l
      
      // get the prior std. dev. for the random coefficients
      vector<Type> PriorSD_bias = exp(log_variance_component_bias);
      // zero off the loop control variable
      LCV = 0;
      for (int l = 0; l < bf_per_res.size(); ++l) {                       // loop through basis function resolutions
        for (int k = 0; k < bf_per_res(l); ++k) {                         // loop through random coefficients within resolution l
          LL_random += dnorm(random(k + LCV), mu, PriorSD_bias(l), true); // add the normal pdf contribution to the log-likelihood
        }
        LCV = LCV + bf_per_res(l);                                        // adjust the loop-control-variable so that next set of K_l coefficients are used
      }
      vector<Type> PriorVar_bias = PriorSD_bias * PriorSD_bias;  // calculate the variance to match that of the VA case
      // Add variance parameter to auto-diff reporting
      ADREPORT(PriorVar_bias);
      
      break;
    }
    case new_latent:
    {
      // additional data specific to IDM accounting for bias with a new secondary latent field
      DATA_IVECTOR(bias_bf_per_res);                  // Integer vector describing the number of basis functions within each spatial resolution
      // Parameters specific to IDM accounting for bias with a secondary latent field
      PARAMETER_VECTOR(log_variance_component_bias);  // log(Standard Dev.) of random coefficients within spatial res., l
      
      // get the prior std. dev. for the random coefficients
      vector<Type> PriorSD_bias = exp(log_variance_component_bias);
      // zero off the loop control variable
      LCV = 0;
      for (int l = 0; l < bias_bf_per_res.size(); ++l) {                  // loop through basis function resolutions
        for (int k = 0; k < bias_bf_per_res(l); ++k) {                    // loop through random coefficients within resolution l
          LL_random += dnorm(random(k + LCV), mu, PriorSD_bias(l), true); // add the normal pdf contribution to the log-likelihood
        }
        LCV = LCV + bias_bf_per_res(l);                                   // adjust the loop-control-variable so that next set of K_l coefficients are used
      }
      vector<Type> PriorVar_bias = PriorSD_bias * PriorSD_bias;  // calculate the variance to match that of the VA case
      // Add variance parameter to auto-diff reporting
      ADREPORT(PriorVar_bias);
      break;
    }
    } // end of basis_type switch
    
    break;
  }
  default:
    error("Model Type not recognised");
  } // end model type switch
