arma::vec movement_function(const int AGE,                  // the current age
                            const int MONTH,                // the current month
                            const int max_cell,             // the max number of cells
                            arma::mat adult_movement_prob,  // adult movement probabilities (from each cell to each cell)
                            arma::cube yearly_pop          // numbers of fish in each cell x month x age
) {
  
  arma::mat pop_before = yearly_pop.slice(AGE).col(MONTH);
  
  arma::vec pop_after = adult_movement_prob.t() * pop_before;
  // adult_movement is the original dataframe, giving the probability of moving *from cell_x to cell_y
  // adult_movement.t() is the transposed version, and is faster to compute with.  

  return pop_after;

}
