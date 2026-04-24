Rcpp::List recruitment_function(
    const int MONTH,            // current spawning month
    const int max_cell,          // number of cells in the grid
    const int max_age,           // oldest age of the species
    const double BHa,            // Beverton-Holt parameter 
    const double BHb,            // Beverton-Holt parameter
    const double PF,             // proportion of females 
    const double ha_scaling,     // hyperallometry scaling factor
    arma::mat maturity,          // age x month matrix of the proportion of the population to be mature over 12 months.
    arma::mat weight,            // age x month matrix of the weight of a fish of each age over months over 12 months.
    arma::vec settlement,   // vector of size max_cell giving probability of recruiting in each cell (based on habitat etc.)
    arma::cube yearly_pop  // numbers of fish in each cell x month x age
) {
  
  // obtain cell x age population for spawning month
  arma::mat pop_month(max_cell, max_age);

  for (int CELL = 0; CELL < max_cell; CELL++)
    for (int AGE = 0; AGE < max_age; AGE++)
      pop_month(CELL, AGE) = yearly_pop(CELL, MONTH, AGE);

  // calculate the spawning biomass for this current month
  arma::vec mature_spwn_month = maturity.col(MONTH);   // % of mature fish of each age at current month
  arma::vec weight_spwn_month = weight.col(MONTH);   // weight of fish of each age at current month

  arma::mat SB_mat(max_cell, max_age);

  for (int AGE = 0; AGE < max_age; AGE++) {
    double age_weight = weight_spwn_month(AGE);
    double ha_weight  = pow(age_weight, ha_scaling); // hyperallometric fecundity per fish
    SB_mat.col(AGE)   = pop_month.col(AGE) * PF * mature_spwn_month(AGE) * ha_weight;
  }
  
  // Sum across ages and cells to get total effective spawning output
  double total_female_SB = arma::accu(SB_mat);
  
  // Standard BH on the hyperallometry-adjusted spawning output
  double tot_recs_before_var = total_female_SB / (BHa + BHb * total_female_SB);
  double tot_recs = tot_recs_before_var * exp(R::rnorm(0, 0.5) - (0.5 * 0.6 * 0.6));
  
  arma::vec settle_recs = settlement * tot_recs;

  return Rcpp::List::create(Rcpp::Named("settle_recs") = settle_recs // vec of each cell's recruits for this spawning month
  );
}
