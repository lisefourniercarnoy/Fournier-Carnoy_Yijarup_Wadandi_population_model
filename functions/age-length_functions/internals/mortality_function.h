Rcpp::List mortality_function(
    // to select current point in time - the master function fills this.
    const int AGE, 
    const int MONTH,
    const int YEAR,
    const int max_cell,
    const int n_lengths,

    // enter parameters
    const double natural_mortality,
    arma::vec weight, // vec of size n_lengths
    arma::mat selectivity, // selectivity-retention of fish in every n_lengths x max_year -- currently a single selectivity for all fleets
    arma::cube current_pop, // cell x length x age
    Rcpp::List fleet_fishing_effort, // from the effort_function - list of size n_fleets, made of vectors of size n_cell
    Rcpp::List fleet_info // a list with objects of different shapes that contain information relevant to the distribution of each fleet.

) {
  
  // constants for this age-month-year
  arma::vec sel           = selectivity.col(YEAR); // this assumes that there is a single selectivity for all fleets
  const double natural_m  = natural_mortality / 12;
  arma::vec wt            = weight;
  arma::mat N             = current_pop.slice(AGE); // max_cell x length of the current age
  int n_fleets            = fleet_fishing_effort.size(); // number of fleets from effort list
  
  
  // 1. calculate the fishing mortality of each fleet
  arma::mat total_f(max_cell, n_lengths, arma::fill::zeros); // cell x lengths
  arma::cube fleet_f(max_cell, n_lengths, n_fleets); // cell x lengths
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    
    arma::vec effort_f    = Rcpp::as<arma::vec>(fleet_fishing_effort[FLEET]);
    
    Rcpp::List current_fleet = fleet_info[FLEET]; // extract information for the current fleet
    arma::vec catchability_now  = Rcpp::as<arma::cube>(current_fleet["catchability"]).slice(YEAR).col(MONTH);
    
    arma::mat F           = (effort_f % catchability_now) * sel.t(); // instantaneous fishing mortality(AGE) = fishing effort(cells) * catchability(cells) * selectivity(all lengths)
    fleet_f.slice(FLEET)  = F;
    total_f += F;  // accumulate total F, as in https://academic.oup.com/icesjms/article/78/6/2043/6317566 -- sum of all fleets' fishing mortalities.
    
  }

  // 2. calculate total survived fish (not fleet-specific)
  arma::mat Z             = total_f + natural_m;
  Z.replace(0, 1e-10); // prevent division by zero
  arma::mat tot_survived  = N % arma::exp(-Z); 
  
  
  // 3. calculate fleet-specific Baranov catch
  Rcpp::List catch_weight_by_fleet(n_fleets); // this is fully age-length structured because we want to know what sizes we're catching
  Rcpp::List catch_numbers_by_fleet(n_fleets); // this is summed across all ages and lengths, to distribute effort next month.
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    // each fleet's share of catch via Baranov, using total Z in denominator
    arma::mat fleet_catch           = N % (fleet_f.slice(FLEET) / Z) % (1 - arma::exp(-Z)); // Baranov catch equation, results in a cell x lengths
    arma::vec fleet_catch_weight    = fleet_catch * wt;
    catch_numbers_by_fleet[FLEET]   = fleet_catch;
    catch_weight_by_fleet[FLEET]    = fleet_catch_weight;
  }

  return Rcpp::List::create(
    Rcpp::Named("fishing_mortality")      = total_f, // cell x length. this is just for checking that fishing mortality is reasonable
    Rcpp::Named("tot_survived")           = tot_survived, // matrix of size max_cell x lengths
    Rcpp::Named("catch_numbers_by_fleet") = catch_numbers_by_fleet, // list of size n_fleets, with matrices of size max_cell x lengths
    Rcpp::Named("catch_weight_by_fleet")  = catch_weight_by_fleet // list of size n_fleets, with vecs of size max_cell
  );
}
