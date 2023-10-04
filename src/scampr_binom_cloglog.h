  // Binomial Generalised Linear Model with complementary log-log link function //
  // Either with or without spatial random effects (the former using a Laplace approx.)

  // Presence-Absence Data
  DATA_MATRIX(X_PA);                        // Fixed effects n_{PA} x p matrix - at the sites
  DATA_VECTOR(Y);                           // Binary vector describing presence/absence - length n_{PA}
  DATA_VECTOR(OFFSET);                      // Fixed effect vector of offset values - length n_{PA}

  // Utility
  DATA_INTEGER(pa_offset);                  // Integer = 0 for no offset term; = 1 for an offset term

  // Initialise the various log-likelihood components
  Type LL_PA = 0.0;

  // Shared by all model types
  vector<Type> Xfixed_PA = X_PA * fixed; // fixed effects at the survey sites (including intercept depending on formula)

  // add in the offset term if present
  if (pa_offset == 1) {
    Xfixed_PA = Xfixed_PA + OFFSET;
  }

  // Initialise the n_{PA} vector of probabiities of absence at the survey sites
  vector<Type> prob_abs(Xfixed_PA.size());

  switch (approx_type){
  case not_sre:
  {
    // probability of absence at the sites  //
    prob_abs = exp(-exp(Xfixed_PA));
    break;
  }
  case variational: // NO CLOSED FORM SOLUTION!
  {
    error("Model Type not recognised");
  }
  case laplace:
  {

    // Initialise the vector resulting from the switch
    vector<Type> Zrandom_PA; // latent field approx. at the survey sites

    // Calculate the random components based on switch
    switch (bf_type){
    case sparse:
    {
      DATA_SPARSE_MATRIX(Z_PA);   // Basis functions n_{PA} x k matrix - at the sites
      Zrandom_PA = Z_PA * random;
      break;
    }
    case dense:
    {
      DATA_MATRIX(Z_PA);          // Basis functions n_{PA} x k matrix - at the sites
      Zrandom_PA = Z_PA * random;
      break;
    }  default:
      error("Basis function matrix type not recognised");
    } // end bf_type switch

    // probability of absence at the sites
    prob_abs = exp(-exp(Xfixed_PA + Zrandom_PA));

    // NOTE: Random coefficient likelihood component is included in scampr.cpp

    break;
  }
  default:
    error("Model Type not recognised");
  } // end model type switch

  // log-likelihood component at the presence/absence sites //

  Type Size = 1;                                        // initialise the number of trials (requires loop due to varying probabilities)
  for (int i = 0; i < X_PA.rows(); i++) {               // loop through survey sites
    LL_PA += dbinom(Y(i), Size, 1 - prob_abs(i), true); // add the binomial pdf contribution to the log-likelihood
  }
