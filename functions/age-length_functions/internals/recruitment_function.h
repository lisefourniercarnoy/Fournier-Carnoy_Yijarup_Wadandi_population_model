Rcpp::List recruitment_function(
    const int MONTH,             // current spawning month
    const int max_cell,          // number of cells in the grid
    const int max_age,           // oldest age of the species
    const double BHa,            // Beverton-Holt parameter 
    const double BHb,            // Beverton-Holt parameter
    const double PF,             // proportion of females 
    const double ha_scaling,     // hyperallometry scaling factor
    const arma::vec maturity,          // vec of size ncell of the proportion of the each length to be mature.
    const arma::vec weight,            // vec of size ncell of the weight of each length.
    const arma::vec settlement,        // vector of size max_cell giving probability of recruiting in each cell (based on habitat etc.)
    const arma::cube current_pop       // numbers of fish in each cell x length x age
) {

  arma::mat SB_age(max_cell, max_age, arma::fill::zeros);
  arma::vec ha_weight = arma::pow(weight, ha_scaling); // vec of size n_lengths — fecundity per fish, adjusted for hyperallometry
  arma::vec maturity_x_haweight = maturity % ha_weight;  // n_lengths vec — spawning weight (not specific to the current pop)
  
  double total_female_SB = 0.0;
  
  for (int AGE = 0; AGE < max_age; AGE++) {
    arma::vec weight_spawning = current_pop.slice(AGE) * maturity_x_haweight; // spawning weight specific to the current pop
    total_female_SB += PF * arma::sum(weight_spawning); // Sum across ages and cells to get total effective spawning output
  }
  
  // Standard BH on the hyperallometry-adjusted spawning output
  double tot_recs_before_var = total_female_SB / (BHa + BHb * total_female_SB);
  
  // we want to add variability in the recruitment. we want the average recruitment to be what the Beverton-Holt equation predicts (variability averaging 1), but some years above, some years below.
  // we exponentiate to make some years *really good* and some years *really bad* for recruitment. exponentiating makes the normal distribution asymmetrical though! exp(1) = 2.72, but exp(-1) = 0.37 (== the average is greater than 1)
  // therefore we do some math to make sure the variability is exponential, but always averaging 1: exp(variability - (sigma^2 / 2))
  double sigma = 0.5; 
  double tot_recs_after_var = tot_recs_before_var * exp(R::rnorm(0, sigma) - ((sigma*sigma)/2)); // add some recruitment variability
  
  arma::vec settle_recs = settlement * tot_recs_after_var;

  return Rcpp::List::create(
    Rcpp::Named("settle_recs") = settle_recs // vec of each cell's recruits for this spawning month
  );
}
