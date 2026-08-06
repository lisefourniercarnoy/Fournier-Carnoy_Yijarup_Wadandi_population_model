Rcpp::List pseudo_effort_function(
    const int MONTH,
    const int max_cell,
    const int n_lengths,
    const int max_age,
    const double natural_mortality,
    const arma::vec& sel_this_year, // vec of size n_lengths
    const arma::cube& current_pop, // cell x length x age
    const std::vector<arma::mat>& catchability_by_fleet // vector of size n_fleet, with matrices of size ncell x months
) {
  
  // 0. obtain constants and full output objects
  
  const double natural_m = natural_mortality / 12; // monthly natural mortality
  int n_fleets = catchability_by_fleet.size();
  
  arma::mat total_f(max_cell, n_lengths, arma::fill::zeros);
  arma::cube fleet_f(max_cell, n_lengths, n_fleets);
  
  
  // 1. calculate pseudo-effort for each fleet
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    arma::vec catchability_now = 100 * catchability_by_fleet[FLEET].col(MONTH);
    arma::mat F = catchability_now * sel_this_year.t(); // here we're not using effort values, because this is a hypothetical fishing mortality
    fleet_f.slice(FLEET) = F;
    total_f += F;
  }
  
  
  // 2. calculate the expected catch from the hypothetical effort
  
  arma::mat Z = total_f + natural_m;
  Z.replace(0, 1e-10); // avoid calculation misbehaviour
 
  Rcpp::List pseudo_catch_by_fleet(n_fleets);
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
  
    arma::vec fleet_catch_total(max_cell, arma::fill::zeros); // ADDED - accumulator across ages
  
    for (int AGE = 0; AGE < max_age; AGE++) {
      arma::mat current_pop_AGE = current_pop.slice(AGE); // cell x length, this age's slice
    
      arma::vec fleet_catch = arma::sum(
        current_pop_AGE % (fleet_f.slice(FLEET) / Z) % (1 - arma::exp(-Z)),
        1  // dim = 1 sums across columns (lengths), returning a max_cell x 1 column vector
      );
      fleet_catch_total += fleet_catch; // sum expected catch over all ages
    }
    
    
    // 3. add a cap to the expected catch
    
    fleet_catch_total = arma::clamp(fleet_catch_total, 0, 10); // this prevents values that are too big from overpowering the distribution of effort.
    
    
    // 4. add uncertainty around the expected catch
    
    double catch_uncertainty = std::pow(R::rnorm(1, 0.3), 2);
    pseudo_catch_by_fleet[FLEET] = fleet_catch_total * std::sqrt(catch_uncertainty); // basically using the absolute value of uncertainty to avoid negative expected catch.
  }
  
  return Rcpp::List::create(
    Rcpp::Named("pseudo_catch") = pseudo_catch_by_fleet
  );
}