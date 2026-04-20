Rcpp::List distribute_effort_function(
    // to select current point in time - the master function fills this.
    const int MONTH, 
    const int YEAR, 
    
    // parameters for the expected catch function
    const int max_cell,
    Rcpp::List expected_catch, 
    Rcpp::List expected_catch_sq,
    
    // fleets
    Rcpp::CharacterVector fleet_names, // c("commercial", "boat_rec", "shore_rec")
    Rcpp::List fleet_info // a list with objects of different shapes that contain information relevant to the distribution of each fleet.

// fleet info should have the following objects:

// utility (cell x access_point x years),
// fishability (cell x month x years), 
// fishing_days (month x access_point x year)
// cell_area_m2 (vector of length n_cell)
// cell_coef (tibble of 4 cols, 1 row)

) {
  
  int n_fleets = fleet_names.size(); // count the number of fleets the calculate
  arma::mat effort_matrix(max_cell, n_fleets); // effort matrix (cells x fleets)
  Rcpp::List fishing_effort_all_fleets(n_fleets); // vector of size n_cells for every fleet, list of size n_fleets
  
  Rcpp::List final_expected_catch(n_fleets); // fleet list of ncell vectors
      
  for (int i = 0; i < n_fleets; i++) { // for each fleet...

    Rcpp::List current_fleet = fleet_info[i]; // extract information for the current fleet
    
    // 1. calculate expected catch 
    arma::vec catchability_now  = Rcpp::as<arma::cube>(current_fleet["catchability"]).slice(YEAR).col(MONTH);
    arma::vec expected_catch_2    = Rcpp::as<arma::vec>(expected_catch[i]) % catchability_now;  // extract ONCE here
    arma::vec expected_catch_sq_2 = arma::pow(expected_catch_2, 2);                             // derive sq from ec
    final_expected_catch[i]     = expected_catch_2;

    // 2. obtain cell utility (different fleets might use cells differently so has to be calculated exogenously
    arma::cube cell_utility       = Rcpp::as<arma::cube>(current_fleet["utility"]); // (cell x access_point x years)
    arma::mat cell_utility_now    = cell_utility.slice(YEAR); // (cell x access_point)
    
    arma::vec cell_area           = Rcpp::as<arma::vec>(current_fleet["cell_area_m2"]); // (vec of length n_cell)

    
    // calculate the coefficients for every access_point
    const int n_access_points = cell_utility_now.n_cols;
    Rcpp::DataFrame coef_values   = Rcpp::as<Rcpp::DataFrame>(current_fleet["coef_values"]); // coefficients from matt's paper
    
    arma::mat utility_calc(max_cell, n_access_points); // announce new matrix (cell x access_point)
    for(int access_point = 0; access_point < n_access_points; access_point++) {
      
      // 2.A. figure out coefficients
      arma::vec cell_coefficent_here_now = // Charlotte's CellCoef
        cell_utility_now.col(access_point) * Rcpp::as<Rcpp::NumericVector>(coef_values["utility"])[0] +
        expected_catch_2                   * Rcpp::as<Rcpp::NumericVector>(coef_values["expected_catch"])[0] +
        expected_catch_sq_2                * Rcpp::as<Rcpp::NumericVector>(coef_values["expected_catch_sq"])[0] + 
        cell_area                          * Rcpp::as<Rcpp::NumericVector>(coef_values["cell_area"])[0];
      
      // arma::vec cell_utility_here_now = arma::exp(cell_coefficent_here_now);
      
      // 2.B. calculate cell utility
      double max_coef = cell_coefficent_here_now.max(); // centering the coefficients to prevent division by zero ??
      arma::vec cell_utility_here_now = arma::exp(cell_coefficent_here_now - max_coef);
      double utility_sum = arma::sum(cell_utility_here_now);
      utility_calc.col(access_point) = (utility_sum > 0) 
        ? cell_utility_here_now / utility_sum 
      : arma::vec(max_cell, arma::fill::value(1.0 / max_cell));
      
      //double utility_sum = arma::sum(cell_utility_here_now);
      
      //utility_calc.col(access_point) = cell_utility_here_now / utility_sum; // calculate each cell's utility from the perspective of the current access point.
    }
    
    // 3. distribute effort
    arma::cube fishing_days         = Rcpp::as<arma::cube>(current_fleet["fishing_days"]); // (month x access_point x years)
    arma::rowvec fishing_days_now   = fishing_days.slice(YEAR).row(MONTH); // (vector of length access_point)
    
    arma::mat fishing_effort        = utility_calc.each_row() % fishing_days_now; // multiply a matrix rows (access_point) with vector (access_points), across all cells
    arma::vec fishing_effort_sum    = arma::sum(fishing_effort, 1); // cells' fishing effort right now, coming from all access points
    fishing_effort_all_fleets[i]    = fishing_effort_sum; // record the fishing effort vector for this fleet
    
    effort_matrix.col(i)            = fishing_effort_sum; // stack all fleet efforts together
  } 
  
  arma::vec total_fishing_effort = arma::sum(effort_matrix, 1); // sum across columns (fleets) to get a vector of size n_cells
  
  return Rcpp::List::create(
    Rcpp::Named("fleet_fishing_effort") = fishing_effort_all_fleets, // list of vectors of size n_cells for every fleet
    Rcpp::Named("total_fishing_effort") = total_fishing_effort // vector of size n_cells for all fleet
  );
} 
