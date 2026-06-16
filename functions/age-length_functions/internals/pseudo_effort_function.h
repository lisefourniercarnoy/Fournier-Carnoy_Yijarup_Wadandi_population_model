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
    // we want to add uncertainty in the expected catch. we want the average recruitment to be what last month's catchwas (uncertainty averaging 1), but some months above, some months below.
    // we exponentiate to make some months *really good* and some months *really bad* for expected catch. exponentiating makes the normal distribution asymmetrical though! exp(1) = 2.72, but exp(-1) = 0.37 (== the average is greater than 1)
    // therefore we do some math to make sure the variability is exponential, but always averaging 1: exp(uncertainty - (sigma^2 / 2))
    double sigma = 0.5;
    double uncertainty = exp(R::rnorm(0, sigma) - ((sigma * sigma) / 2)); // exp() because
    arma::vec temp = past_month_cpue * uncertainty; // transform into *expected catch* by adding uncertainty - store in a temporary object
    expected_catch[i] = temp; // transform into *expected catch* by adding uncertainty
    expected_catch_sq[i] = arma::pow(temp, 2); // square it

  }

 return Rcpp::List::create(Rcpp::Named("expected_catch") = expected_catch,
                           Rcpp::Named("expected_catch_sq") = expected_catch_sq);
 }
