#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

// THIS FUNCTION IS THE SAME AS CHARLOTTE'S BUT WAS ANNOTATED BY LISE (for her brain)
// RECRUITMENT FUNCTION WAS CHANGED (referred to the spawning month of fish, changed from 9 -october- to 10 -november-)

// [[Rcpp::export]]
arma::vec movementfunc_cpp(const int AGE, 
                           const int MONTH, 
                           const int MaxCell, 
                           arma::mat AdultMove, 
                           arma::cube YearlyTotal) {
  
  // declare object types
  int Cell_cl;
  int Cell_rw;
  double Pop;
    
  arma::vec All_Movers(MaxCell);
  arma::mat Pop2(MaxCell,MaxCell);
  
  for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
    for (Cell_cl=0; Cell_cl<MaxCell; Cell_cl++) { 
      
      Pop = YearlyTotal(Cell_rw,MONTH,AGE);
      
      Pop2(Cell_rw,Cell_cl) = Pop2(Cell_rw,Cell_cl) + (AdultMove(Cell_rw,Cell_cl) * Pop); 
    } 
  } // this loop calculates how many fish will move from each cell to each other cell (done for every age and every month)
  
  for (Cell_cl=0; Cell_cl<MaxCell; Cell_cl++) { 
    for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
      All_Movers(Cell_cl)=All_Movers(Cell_cl) + Pop2(Cell_rw,Cell_cl);
    }
  } // this loop sums how many fish are in each cell resulting from movement from all other cells (for each age, and each month)
    
  return All_Movers;
}

// [[Rcpp::export]]
Rcpp::List mortalityfunc_cpp(const int AGE, 
                             const int MaxCell, 
                             const int MONTH, 
                             const int YEAR, 
                             const double NatMort,
                             arma::mat Weight, 
                             arma::cube Selectivity, 
                             arma::cube YearlyTotal, 
                             arma::cube Effort) {
  // declare object types
  int Cell_rw;
  arma::vec tot_survived(MaxCell);
  arma::vec tot_died(MaxCell);
  arma::vec caught(MaxCell);
  arma::vec finite_f(MaxCell);
  arma::vec caught_weight(MaxCell);
  arma::vec Z(MaxCell);
  arma::vec YPR(MaxCell);
  arma::vec Baranov(MaxCell);

  for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
    
    tot_survived(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) * ((exp(-Effort(Cell_rw,MONTH,YEAR) * Selectivity(AGE,MONTH,YEAR))) * exp(-NatMort/12));
    caught(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) * (1-(exp(-Effort(Cell_rw,MONTH,YEAR) * Selectivity(AGE,MONTH,YEAR))));
    tot_died(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) - tot_survived(Cell_rw);
    finite_f(Cell_rw) = Effort(Cell_rw, MONTH, YEAR) * Selectivity(AGE,MONTH,YEAR); // This isn't actually finite but I can't be bothered changing its name
    
    Z(Cell_rw) = (Effort(Cell_rw, MONTH,YEAR) * Selectivity(AGE,MONTH,YEAR)) + NatMort/12;
    YPR(Cell_rw) = tot_survived(Cell_rw)*((Effort(Cell_rw,MONTH,YEAR) * Selectivity(AGE,MONTH,YEAR))/Z(Cell_rw)) * (1-(exp(-Z(Cell_rw)))) * Weight(AGE, MONTH);
    
    Baranov(Cell_rw) = tot_survived(Cell_rw) * (finite_f(Cell_rw)/Z(Cell_rw)) * (1-(exp(-Z(Cell_rw))));
    caught_weight(Cell_rw) = Baranov(Cell_rw) * Weight(AGE, MONTH);
  } // this loop calculates mortality for every age group and every month and every cell

  return Rcpp::List::create(Rcpp::Named("tot_survived") = tot_survived,
                            Rcpp::Named("tot_died") = tot_died,
                            Rcpp::Named("Catch") = caught,
                            Rcpp::Named("finite_f") = finite_f,
                            Rcpp::Named("Catch_weight") = caught_weight,
                            Rcpp::Named("YPR") = YPR,
                            Rcpp::Named("Age_Baranov") = Baranov);
}


// [[Rcpp::export]]
Rcpp::List recruitmentfunc_cpp(const int MaxCell, 
                               const int MaxAge, 
                               const double BHa, 
                               const double BHb, 
                               const double PF, 
                               arma::mat Mature, 
                               arma::mat Weight, 
                               arma::vec Settlement, 
                               arma::cube YearlyTotal) {
  // declare object types
  int AGE;
  int Cell_rw;
  double Fem_adults; 
  double tot_recs;
  arma::vec SB(MaxCell);
  arma::vec recs(MaxAge);
  arma::vec recs_variable(MaxAge);
  arma::vec settle_recs(MaxCell);
  arma::vec TotFemSB(MaxAge);
  
  for (AGE=0; AGE<MaxAge; AGE++) { 
    for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
      Fem_adults = YearlyTotal(Cell_rw,10,AGE) * PF; // month = November Female adults of that age class
      SB(Cell_rw) = Fem_adults * (Mature(AGE,10)) * (Weight(AGE,10)); // Gives us biomass
    }
    TotFemSB(AGE) = sum(SB); // adding up across the cells - all biomass for one age group
    
    recs(AGE) = (TotFemSB(AGE) / (BHa + BHb*TotFemSB(AGE))); // calculate the number of recs from that age group
    recs_variable(AGE) = recs(AGE); 
  } // this loop calculates how many recruits there are based on the female population in October

  tot_recs = (sum(recs_variable)) * exp(R::rnorm(0, 0.6)-(0.5*0.6*0.6));

  for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
    settle_recs(Cell_rw) = Settlement(Cell_rw) * tot_recs; 
  } // this loop calculates how many fish are recruited from settlement
  
  return Rcpp::List::create(Rcpp::Named("settle_recs") = settle_recs,
                            Rcpp::Named("BH_recs") = recs,
                            Rcpp::Named("Fem_adults") = Fem_adults, // All these go once fixed
                            Rcpp::Named("SB") = SB,
                            Rcpp::Named("TotFemSB") = TotFemSB,
                            Rcpp::Named("tot_recs") = tot_recs);
  
  }

// [[Rcpp::export]]
Rcpp::List RunModelfunc_cpp(const int YEAR, 
                            const int MaxAge, 
                            const int MaxYear, 
                            const int MaxCell, 
                            const double NatMort, 
                            const double BHa, 
                            const double BHb, 
                            const double PF, 
                            arma::mat AdultMove, 
                            arma::mat Mature, 
                            arma::mat Weight, 
                            arma::vec Settlement,
                            arma::cube YearlyTotal, 
                            arma::cube Select, 
                            arma::cube Effort) {
  
  // declare object types
  int MONTH;
  int AGE;
  int Cell_rw;
  Rcpp::List Res;
  Rcpp::List recruits;
  Rcpp::NumericVector Fish_catch(MaxCell);
  Rcpp::NumericVector tot_mort(MaxCell);
  Rcpp::NumericVector caught(MaxCell);
  Rcpp::NumericVector caught_weight(MaxCell);
  Rcpp::NumericVector tot_survived(MaxCell);
  Rcpp::NumericVector finite_f(MaxCell);
  Rcpp::NumericVector tot_died(MaxCell);
  Rcpp::NumericVector Weight_catch(MaxCell);
  Rcpp::NumericVector All_Movers(MaxCell);
  Rcpp::NumericVector settle_recs(MaxCell);
  Rcpp::NumericVector BH_recs(MaxCell);
  Rcpp::NumericVector Fem_SB(MaxAge);
  Rcpp::NumericVector Z(MaxCell);
  Rcpp::NumericVector YPR(MaxCell);
  Rcpp::NumericVector YPR_month(MaxCell);
  
  arma::cube month_catch(MaxCell, 12, MaxAge);
  arma::cube month_catch_weight(MaxCell, 12, MaxAge);
  arma::cube age_catch(MaxCell, 12, MaxAge);
  arma::cube tot_mort_all(MaxCell, 12, MaxAge);
  arma::cube age_dead(MaxCell, 12, MaxAge);
  arma::cube age_survived(MaxCell, 12, MaxAge);
  arma::cube month_YPR(MaxCell, 12, MaxAge);
  
  for (MONTH=0; MONTH<12; MONTH++) { 
    // the function is ordered as follows: 1. fish move to another cell, 2. fish die, 3. fish reproduce
    
    // 1. movement
    
    for (AGE=0; AGE<MaxAge; AGE++) { 
      
      All_Movers = movementfunc_cpp(AGE, MONTH, MaxCell, AdultMove, YearlyTotal);

      for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {

        YearlyTotal(Cell_rw,MONTH,AGE) = All_Movers(Cell_rw);
      }
      
    } // this loop moves fish, for every cell and across all ages
    
    // 2. mortality
    
    for (AGE=0; AGE<MaxAge; AGE++) {
      
      if (MONTH==11) {
        if (AGE<(MaxAge-1)) {
          
          Res = mortalityfunc_cpp(AGE, MaxCell, MONTH, YEAR, NatMort, Weight, Select, YearlyTotal, Effort);
          
          tot_survived = Res["tot_survived"];
          tot_died = Res["tot_died"];
          caught = Res["Age_Baranov"];
          caught_weight = Res["Catch_weight"];
          YPR_month = Res["YPR"];
          
          
          for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
            YearlyTotal(Cell_rw,0,AGE+1) = tot_survived(Cell_rw);
          }
        } // this loop graduates fish up to the next age class if they've survived to December
        
      } else if (MONTH<11){ // Month==11
        
        Res = mortalityfunc_cpp(AGE, MaxCell, MONTH, YEAR, NatMort, Weight, Select, YearlyTotal, Effort);
        
        tot_survived = Res["tot_survived"];
        tot_died = Res["tot_died"];
        caught = Res["Age_Baranov"];
        caught_weight = Res["Catch_weight"];
        YPR_month = Res["YPR"];
        
        for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
          YearlyTotal(Cell_rw,MONTH+1,AGE) = tot_survived(Cell_rw);
        }
      } // this loop applies the mortality function to fish in each age class and each month January-November, and ages up December survivors
      
      // setting other outputs 
      for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
        month_catch(Cell_rw,MONTH,AGE) = caught(Cell_rw);
        month_catch_weight(Cell_rw,MONTH,AGE) = caught_weight(Cell_rw);
        age_dead(Cell_rw,MONTH,AGE) = tot_died(Cell_rw);
        age_survived(Cell_rw,MONTH,AGE) = tot_survived(Cell_rw);
        month_YPR(Cell_rw, MONTH, AGE) = YPR_month(Cell_rw); // yield per recruit
      }
    } // this loop calculates extra metrics (inc. fishing and other mortality metrics)
    
    // 3. Recruitment
    
    if (MONTH==9) { // October = spawning
      recruits = recruitmentfunc_cpp(MaxCell, MaxAge, BHa, BHb, PF, Mature, Weight, Settlement, YearlyTotal);
      
      settle_recs = recruits["settle_recs"];
      BH_recs = recruits["BH_recs"];
      Fem_SB = recruits["TotFemSB"];

      for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
        YearlyTotal(Cell_rw,0,0) = settle_recs(Cell_rw);
      }
    } // this loop applies the recruitment function to the population in the spawning month
  }
  
  return Rcpp::List::create(Rcpp::Named("YearlyTotal") = YearlyTotal,
                            Rcpp::Named("settle_recs") = settle_recs,
                            Rcpp::Named("BH_recs") = BH_recs,
                            Rcpp::Named("month_catch") = month_catch,
                            Rcpp::Named("month_catch_weight") = month_catch_weight,
                            Rcpp::Named("age_died") = age_dead,
                            Rcpp::Named("age_survived") = age_survived,
                            Rcpp::Named("Monthly_YPR_age") = month_YPR,
                            Rcpp::Named("Fem_SB") = Fem_SB);
  
} // End function