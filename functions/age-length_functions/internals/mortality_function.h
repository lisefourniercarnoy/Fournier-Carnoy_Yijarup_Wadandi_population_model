Rcpp::List mortality_function(
    // to select current point in time - the master function fills this.
    const int AGE, 
    const int MONTH,
    const int max_cell,
    const int n_lengths,

    // enter parameters
    const double natural_mortality,
    const arma::vec& weight,
    const arma::vec& sel_this_year, // selectivity-retention of fish in this year (extracted in the master function). currently a single selectivity for all fleets
    const arma::mat& current_pop_AGE, // cell x length
    const Rcpp::List& fleet_fishing_effort, // from the effort_function - list of size n_fleets, made of vectors of size n_cell
    const std::vector<arma::mat>& catchability_by_fleet // extracted before the MONTH loop in the master function
) {
  
  // constants for this age-month-year
  const double natural_m = natural_mortality / 12;
  int n_fleets = fleet_fishing_effort.size(); // number of fleets from effort list
 
  // 1. calculate the fishing mortality of each fleet
  arma::mat total_f(max_cell, n_lengths, arma::fill::zeros); // cell x lengths
  arma::cube fleet_f(max_cell, n_lengths, n_fleets); // cell x lengths
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {

    arma::vec effort_f = Rcpp::as<arma::vec>(fleet_fishing_effort[FLEET]);
    arma::vec catchability_now = catchability_by_fleet[FLEET].col(MONTH);
    arma::mat F = (effort_f % catchability_now) * sel_this_year.t(); // instantaneous fishing mortality(AGE) = fishing effort(cells) * catchability(cells) * selectivity(all lengths)
    fleet_f.slice(FLEET) = F;
    total_f += F;  // accumulate total F, as in https://academic.oup.com/icesjms/article/78/6/2043/6317566 -- sum of all fleets' fishing mortalities.
    
  }

  // 2. calculate total survived fish (not fleet-specific)
  arma::mat Z = total_f + natural_m;
  Z.replace(0, 1e-10); // prevent division by zero
  arma::mat tot_survived = current_pop_AGE % arma::exp(-Z); 

  // 3. calculate fleet-specific Baranov catch
  Rcpp::List catch_weight_by_fleet(n_fleets); // this is fully age-length structured because we want to know what sizes we're catching
  Rcpp::List catch_numbers_by_fleet(n_fleets); // this is summed across all ages and lengths, to distribute effort next month.
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    // each fleet's share of catch via Baranov, using total Z in denominator
    arma::mat fleet_catch = current_pop_AGE % (fleet_f.slice(FLEET) / Z) % (1 - arma::exp(-Z)); // Baranov catch equation, results in a cell x lengths
    arma::vec fleet_catch_weight = fleet_catch * weight;
    catch_numbers_by_fleet[FLEET] = fleet_catch;
    catch_weight_by_fleet[FLEET] = fleet_catch_weight;
   
  }

  return Rcpp::List::create(
    Rcpp::Named("fishing_mortality")      = total_f, // cell x length. this is just for checking that fishing mortality is reasonable
    Rcpp::Named("fishing_mortality_by_fleet") = fleet_f,   // cell x length x fleet
    Rcpp::Named("tot_survived")           = tot_survived, // matrix of size max_cell x lengths
    Rcpp::Named("catch_numbers_by_fleet") = catch_numbers_by_fleet, // list of size n_fleets, with matrices of size max_cell x lengths
    Rcpp::Named("catch_weight_by_fleet")  = catch_weight_by_fleet // list of size n_fleets, with vecs of size max_cell
  );
}
