arma::mat movement_function(const int AGE, // the current age, set in the master function loop
                            const int MONTH, // current month, for spawning movement
                            const int max_cell, // the max number of cells
                            const arma::mat& adult_movement_prob, // adult movement probabilities (>375mm) (from each cell to each cell)
                            const arma::mat& juv_movement_prob, // juvenile movement probabilities (<375mm) (from each cell to each cell)
                            const arma::mat& spawn_movement_prob, // adult spawner movement probabilities (>625mm) (from each cell to each cell)
                            const arma::mat& current_pop_AGE // numbers of fish in each cell x length x age
) {
  
  const int adult_length_bin = 9 - 1; // length bin n.8 is the one with a midpoint of 375  (we'll call that the last juvenile bin). zero-index it!
  const int spawning_length_bin = 13 - 1; // length bin n.13 is the one with a midpoint of 625  (we'll call that the first adult-spawner bin). zero-index it!
  
  arma::mat pop_after(current_pop_AGE.n_rows, current_pop_AGE.n_cols, arma::fill::zeros);
  
  // adult movement
  pop_after.cols(adult_length_bin, pop_after.n_cols - 1) = // for all adult length bins ...
    adult_movement_prob.t() * current_pop_AGE.cols(adult_length_bin, current_pop_AGE.n_cols - 1); // ... move them
  
  // juvenile movement
  pop_after.cols(0, adult_length_bin - 1) = // for all juvie length bins ...
    juv_movement_prob.t() * current_pop_AGE.cols(0, adult_length_bin - 1); // ... move them
  
  // spawning movement
  if (MONTH >= 9) { // for months October onwards,
    pop_after.cols(spawning_length_bin, pop_after.n_cols - 1) = // for all adult-spawner bins ...
      spawn_movement_prob.t() * current_pop_AGE.cols(spawning_length_bin, current_pop_AGE.n_cols - 1); // ... move them (this overrides the adult movement just before)
  }
  
  // adult_movement is the original dataframe, giving the probability of moving *from cell_x to cell_y
  // adult_movement.t() is the transposed version, and is faster to compute with.  

  return pop_after;

}
