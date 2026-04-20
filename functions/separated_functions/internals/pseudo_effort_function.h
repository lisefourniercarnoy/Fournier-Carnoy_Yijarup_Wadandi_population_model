Rcpp::List pseudo_effort_function(
    
    Rcpp::List catch_weight_by_fleet, // vectors (ncell) of catch weight, list object per fleet -- made in the mortality function of the previous month
    Rcpp::List fleet_fishing_effort // vectors (ncell) of fishing effort, list object per fleet -- made in the effort function of the previous month
    
    ) {
  
  const int n_fleets = catch_weight_by_fleet.size();
  Rcpp::List expected_catch(n_fleets);
  Rcpp::List expected_catch_sq(n_fleets);
  
  for (int i = 0; i < n_fleets; i++) { // for each fleet...
   
   // obtain previous month's catch, effort, and calculate CPUE
    arma::vec catch_this_fleet = catch_weight_by_fleet[i];
    arma::vec fishing_effort_this_fleet = fleet_fishing_effort[i];
    fishing_effort_this_fleet.replace(0, 1e-10); // prevent division by zero
    
    
    arma::vec past_month_cpue = catch_this_fleet / fishing_effort_this_fleet;
    past_month_cpue.replace(arma::datum::inf, 0); // prevent issues
    past_month_cpue.replace(arma::datum::nan, 0); // prevent issues
    
    // add uncertainty 
    double uncertainty = R::rnorm(0, 2);
    arma::vec temp = past_month_cpue * uncertainty; // transform into *expected catch* by adding uncertainty - store in a temporary object
    expected_catch[i] = temp; // transform into *expected catch* by adding uncertainty
    expected_catch_sq[i] = arma::pow(temp, 2); // square it

  }

 return Rcpp::List::create(Rcpp::Named("expected_catch") = expected_catch,
                           Rcpp::Named("expected_catch_sq") = expected_catch_sq);
 }
