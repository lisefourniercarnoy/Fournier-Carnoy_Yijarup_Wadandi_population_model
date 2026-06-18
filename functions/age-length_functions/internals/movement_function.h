arma::mat movement_function(const int AGE,                  // the current age, set in the master function loop
                            const int max_cell,             // the max number of cells
                            const arma::mat& adult_movement_prob,   // adult movement probabilities (from each cell to each cell)
                            const arma::mat& current_pop_AGE    // numbers of fish in each cell x length x age
) {
  
  arma::mat pop_after = adult_movement_prob.t() * current_pop_AGE;
  // adult_movement is the original dataframe, giving the probability of moving *from cell_x to cell_y
  // adult_movement.t() is the transposed version, and is faster to compute with.  

  return pop_after;

}
