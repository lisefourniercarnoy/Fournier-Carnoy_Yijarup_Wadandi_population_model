
Rcpp::List recruitment_function(
    const int MONTH,            // current spawning month
    const int max_cell,          // number of cells in the grid
    const int max_age,           // oldest age of the species
    const double BHa,            // Beverton-Holt parameter 
    const double BHb,            // Beverton-Holt parameter
    const double PF,             // proportion of females 
    const double ha_scaling,     // hyperallometry scaling factor
    arma::mat maturity,          // age x month matrix of the proportion of the population to be mature over 12 months.
    arma::mat weight,            // age x month matrix of the weight of a fish of each age over months over 12 months.
    arma::vec settlement,   // vector of size max_cell giving probability of recruiting in each cell (based on habitat etc.)
    arma::cube yearly_pop  // numbers of fish in each cell x month x age
) {
  
  // obtain cell x age population for spawning month
  arma::mat pop_month(max_cell, max_age);
  
  for (int CELL = 0; CELL < max_cell; CELL++)
    for (int AGE = 0; AGE < max_age; AGE++)
      pop_month(CELL, AGE) = yearly_pop(CELL, MONTH, AGE);
  
  // calculate the spawning biomass for this current month
  arma::vec mature_spwn_month = maturity.col(MONTH);   // % of mature fish of each age at current month
  arma::vec weight_spwn_month = weight.col(MONTH);   // weight of fish of each age at current month
  
  arma::mat SB_mat(max_cell, max_age);
  for (int AGE = 0; AGE < max_age; AGE++)
    SB_mat.col(AGE) = pop_month.col(AGE) * PF * mature_spwn_month(AGE) * weight_spwn_month(AGE); // spawning biomass of each age
  arma::vec total_female_SB = arma::sum(SB_mat, 1).t(); // total female spawning biomass for all ages.
  
  // calculate Beverton-Holt recruitment per age
  arma::vec recs(max_age);
  for (int AGE = 0; AGE < max_age; AGE++)
    recs(AGE) = pow(total_female_SB(AGE), ha_scaling) / (BHa + BHb * total_female_SB(AGE)); // this calculates the recruits produced by each age group (given that bigger fish produce disproportionately more babies)
  
  arma::vec recs_variable = recs;
  double tot_recs = arma::sum(recs_variable) * exp(R::rnorm(0, 0.6) - (0.5 * 0.6 * 0.6)); // add some variability in the recruitment so it leaves some things up to environmental variability.
  
  // distribute recruits by settlement probability
  arma::vec settle_recs = settlement * tot_recs;
  
  return Rcpp::List::create(Rcpp::Named("settle_recs") = settle_recs // vec of each cell's recruits for this spawning month
                            //Rcpp::Named("BH_recs") = recs // not sure i need this?
  );
}







// charlotte's function below

// // [[Rcpp::export]]
// Rcpp::List recruitmentfunc_cpp(const int MaxCell, 
//                                const int MaxAge, 
//                                const double BHa, 
//                                const double BHb, 
//                                const double PF, 
//                                const double ha_scaling, // hyperallometry scaling value
//                                arma::mat Mature, 
//                                arma::mat Weight, 
//                                arma::vec Settlement,
//                                arma::cube YearlyTotal) {
//   // declare object types
//   int AGE;
//   int Cell_rw;
//   double Fem_adults; 
//   double tot_recs;
//   arma::vec SB(MaxCell);
//   arma::vec recs(MaxAge);
//   arma::vec recs_variable(MaxAge);
//   arma::vec settle_recs(MaxCell);
//   arma::vec TotFemSB(MaxAge);
//   
//   for (AGE=0; AGE<MaxAge; AGE++) { 
//     for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
//       Fem_adults = YearlyTotal(Cell_rw,10,AGE) * PF; // month = November Female adults of that age class
//       SB(Cell_rw) = Fem_adults * (Mature(AGE,10)) * (Weight(AGE,10)); // Gives us biomass
//     }
//     TotFemSB(AGE) = sum(SB); // adding up across the cells - all biomass for one age group
//     
//     recs(AGE) = ((pow(TotFemSB(AGE), ha_scaling)) / (BHa + BHb*TotFemSB(AGE))); // calculate the number of recs from that age group, with hyperallometry
//     recs_variable(AGE) = recs(AGE); 
//   } // this loop calculates how many recruits there are based on the female population in October
//   
//   tot_recs = (sum(recs_variable)) * exp(R::rnorm(0, 0.6)-(0.5*0.6*0.6));
//   
//   for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
//     settle_recs(Cell_rw) = Settlement(Cell_rw) * tot_recs; 
//   } // this loop calculates how many fish are recruited from settlement
//   
//   return Rcpp::List::create(Rcpp::Named("settle_recs") = settle_recs,
//                             Rcpp::Named("BH_recs") = recs,
//                             Rcpp::Named("Fem_adults") = Fem_adults, // All these go once fixed
//                             Rcpp::Named("SB") = SB,
//                             Rcpp::Named("TotFemSB") = TotFemSB,
//                             Rcpp::Named("tot_recs") = tot_recs);
//   
// }