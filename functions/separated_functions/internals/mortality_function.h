Rcpp::List mortality_function(
    // to select current point in time - the master function fills this.
    const int AGE, 
    const int MONTH,
    const int YEAR,
    const int max_cell,

    // enter parameters
    const double natural_mortality,
    arma::mat weight, // age x month
    arma::cube selectivity, // age x month x year
    arma::cube yearly_pop, // cell x month x age
    Rcpp::List fleet_fishing_effort, // from the effort_function - list of size n_fleets, made of vectors of size n_cell
    Rcpp::List fleet_info // a list with objects of different shapes that contain information relevant to the distribution of each fleet.

) {
  
  // constants for this age-month-year
  const double sel        = selectivity(AGE, MONTH, YEAR); // this assumes that there is a single selectivity for all fleets
  const double natural_m  = natural_mortality / 12;
  const double wt         = weight(AGE, MONTH);
  arma::vec N             = yearly_pop.slice(AGE).col(MONTH);   // (cell) vector for this AGE/MONTH
  int n_fleets            = fleet_fishing_effort.size();   // number of fleets from effort list
  
  
  // 1. calculate the fishing mortality of each fleet
  arma::vec total_f(max_cell, arma::fill::zeros);
  arma::mat fleet_f(max_cell, n_fleets);
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    
    arma::vec effort_f    = Rcpp::as<arma::vec>(fleet_fishing_effort[FLEET]);
    
    Rcpp::List current_fleet = fleet_info[FLEET]; // extract information for the current fleet
    arma::vec catchability_now  = Rcpp::as<arma::cube>(current_fleet["catchability"]).slice(YEAR).col(MONTH);
    
    arma::vec F           = effort_f % catchability_now * sel; // instantaneous fishing mortality(AGE) = fishing effort * catchability * selectivity(AGE)
    fleet_f.col(FLEET)    = F;
    total_f += F;  // accumulate total F, as in https://academic.oup.com/icesjms/article/78/6/2043/6317566 -- sum of all fleets' fishing mortalities.
    
  }

  // 2. calculate total survived fish (not fleet-specific)
  arma::vec Z             = total_f + natural_m;
  Z.replace(0, 1e-10); // prevent division by zero
  arma::vec tot_survived  = N % arma::exp(-Z); 
  
  
  // 3. calculate fleet-specific Baranov catch
  Rcpp::List catch_weight_by_fleet(n_fleets);
  Rcpp::List catch_numbers_by_fleet(n_fleets);
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    // each fleet's share of catch via Baranov, using total Z in denominator
    arma::vec fleet_catch           = N % (fleet_f.col(FLEET) / Z) % (1 - arma::exp(-Z)); // Baranov catch equation,
    arma::vec fleet_catch_weight    = fleet_catch * wt;
    catch_numbers_by_fleet[FLEET]   = fleet_catch;
    catch_weight_by_fleet[FLEET]    = fleet_catch_weight;
  }

  return Rcpp::List::create(
    Rcpp::Named("fishing_mortality")      = total_f, // vec(ncell) this is just for checking that fishing mortality is reasonable
    Rcpp::Named("tot_survived")           = tot_survived,
    Rcpp::Named("catch_numbers_by_fleet") = catch_numbers_by_fleet,
    Rcpp::Named("catch_weight_by_fleet")  = catch_weight_by_fleet
  );
}
