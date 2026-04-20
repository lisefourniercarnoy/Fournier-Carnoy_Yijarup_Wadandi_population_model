arma::vec movement_function(const int AGE,                 // the current age
                            const int MONTH,               // the current month
                            const int max_cell,             // the max number of cells
                            arma::mat adult_movement_prob, // adult movement probabilities (from each cell to each cell)
                            arma::cube yearly_pop_total    // numbers of fish in each cell x month x age
) {
  
  arma::vec pop_before(max_cell);
  for (int i = 0; i < max_cell; i++) {
    pop_before(i) = yearly_pop_total(i, MONTH, AGE);
  }
  // pop_before = age 1, month 1, number of fish in all cells
  // Note: check your cube indexing dim order matches (cell, month, age)
  
  arma::vec pop_after = adult_movement_prob.t() * pop_before;
  // adult_movement is the original dataframe, giving the probability of moving *from cell_x to cell_y
  // adult_movement.t() is the transposed version, and is faster to compute with.
  
  return pop_after;
}
