#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]


#include "internals/effort_function.h"
#include "internals/mortality_function.h"
#include "internals/pseudo_effort_function.h"
#include "internals/recruitment_function.h"
#include "internals/movement_function.h"

// [[Rcpp::export]]
Rcpp::List run_full_model_function(
    
    const int YEAR,
    const int max_cell,
    const int max_age,
    const int n_lengths, // number of length classes
    const int max_year,
    
    arma::cube current_pop, // number of fish in each cell x length x age
    arma::vec weight, // vec of size n_lengths
    arma::mat selectivity, // selectivity-retention of fish in every n_lengths x max_year
    arma::mat age_transition, // age-transition matrix, n_lengths x n_lengths
    const double natural_mortality, // constant
    
    arma::vec spawning_months,
    const double BHa,
    const double BHb,
    const double PF,
    const double ha_scaling,
    arma::vec maturity, // vec of size n_lengths
    arma::vec settlement, // vec of size max_cell
    
    arma::mat adult_movement_prob, // max_cell x max_cell
    
    Rcpp::CharacterVector fleet_names,
    Rcpp::List fleet_info
) {
  
  
  // 0a. create master objects ------------------------------------------------
  
  // these are objects that only need to be computed once per year. this makes the function go faster
  int n_fleets = fleet_names.size();
  
  arma::vec   january_recruits(max_cell, arma::fill::zeros); // this is for a sum of recruits of all spawning months. gets renewed every year.
  arma::mat   master_SSB(max_cell, 12, arma::fill::zeros);        // SSB per cell x month (filled only in spawning months)
  arma::mat   master_F(max_cell, n_lengths, arma::fill::zeros); // just to test whether fishing mortality is realistic
  
  arma::mat   master_F_numerator(max_cell, n_lengths, arma::fill::zeros);
  arma::mat   master_F_denominator(max_cell, n_lengths, arma::fill::zeros);
  
  
  arma::cube  next_pop(max_cell, n_lengths, max_age, arma::fill::zeros); // object to put survivors and recruits in 

  // below are objects that aren't used in the functions, but are needed to store outputs and use in R for checks etc.
  Rcpp::List  catch_weight(n_fleets);
  Rcpp::List  fleet_fishing_effort(n_fleets);
  arma::cube              master_effort_by_fleet(max_cell, 12, n_fleets, arma::fill::zeros); // format (cell, month, fleet) 
  arma::field<arma::mat>  master_catch_weight_total(n_fleets); // format cell x age x fleet
  arma::field<arma::cube> master_catch_number_total(n_fleets); // format cell x length x age x fleet
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    master_catch_number_total(FLEET) = arma::cube(max_cell, n_lengths, max_age, arma::fill::zeros);
    master_catch_weight_total(FLEET) = arma::mat(max_cell, max_age, arma::fill::zeros);
    
    catch_weight[FLEET] = arma::vec(max_cell, arma::fill::ones);
    fleet_fishing_effort[FLEET] = arma::vec(max_cell, arma::fill::ones);
  } // fill this with zeroes to start with, so that they're the right dimensions.
  
  arma::field<arma::cube> master_pop_survived(12); // format cell x length x age x month
  for (int m = 0; m < 12; m++) {
    master_pop_survived(m) = arma::cube(max_cell, n_lengths, max_age, arma::fill::zeros);
  } // fill things with zeroes
  

  // 0b. extract objects that need to be calculated only once per year --------

  arma::vec sel_this_year = selectivity.col(YEAR);
  
  std::vector<arma::mat> catchability_by_fleet(n_fleets);
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    Rcpp::List fl = fleet_info[FLEET];
    catchability_by_fleet[FLEET] = Rcpp::as<arma::cube>(fl["catchability"]).slice(YEAR);
    // now a max_cell x 12 matrix; index .col(MONTH) inside the loop
  }
  std::vector<bool> is_spawn_month(12, false);
  for (arma::uword k = 0; k < spawning_months.n_elem; k++)
    is_spawn_month[static_cast<int>(spawning_months(k)) - 1] = true;
    
  
  // START OF THE FUNCTION HERE ===============================================
  
  for (int MONTH = 0; MONTH < 12; MONTH++) { 
    
    Rcpp::Rcout << "Year: " << YEAR+1 << ", Month: " << MONTH+1 << std::endl;
    
    // 1. figure out where high CPUE was in the previous month ------------------
    
    Rcpp::List expected_catch(n_fleets); // (list n_fleet, obj vec ncell)
    Rcpp::List expected_catch_sq(n_fleets); // (list n_fleet, obj vec ncell)
    
    if (YEAR == 0 && MONTH == 0) { // for the first month of the first year, no information about catch, so everywhere is the same.
      for (int FLEET = 0; FLEET < fleet_names.size(); FLEET++) {
        expected_catch[FLEET] = arma::vec(max_cell, arma::fill::ones); // (list n_fleet, obj vec ncell)
        expected_catch_sq[FLEET] = arma::vec(max_cell, arma::fill::ones); // (list n_fleet, obj vec ncell)
        
      } 
    } // for all fleets for year1 month1, fill the expected catch (and expected catch squared) with zeroes, since we don't have any previous month CPUE.
    
    else {
      
      Rcpp::List pseudo = pseudo_effort_function(
        catch_weight, // list n_fleet, obj vec ncell -- made in the mortality function of the previous month
        fleet_fishing_effort // vectors (ncell) of fishing effort, list object per fleet -- made in the effort function of the previous month
      );
      
      expected_catch = pseudo["expected_catch"]; // (list n_fleet, obj vec ncell)
      expected_catch_sq  = pseudo["expected_catch_sq"]; // (list n_fleet, obj vec ncell)
      
    } // for all months and years after y1 m1, use the previous month's catch and effort distribution to obtain pseudo catch for this month
      
      Rcpp::Rcout << "fisher strategy done - ";
    
    
    // 2. distribute effort accordingly ---------------------------------------

    Rcpp::List effort_distribution_output = distribute_effort_function(MONTH,
                                                                       YEAR,
                                                                       max_cell,
                                                                       expected_catch,
                                                                       fleet_fishing_effort,
                                                                       fleet_names,
                                                                       fleet_info);
    
    // outputs needed to make current month's effort distribution visible to next month's (pseudo effort function)
    fleet_fishing_effort = Rcpp::as<Rcpp::List>(effort_distribution_output["fleet_fishing_effort"]); // list n_fleet, obj vec ncell
    arma::vec total_fishing_effort = Rcpp::as<arma::vec>(effort_distribution_output["total_fishing_effort"]); // vec ncell

    // also fill the master outputs
    for (int FLEET = 0; FLEET < n_fleets; FLEET++) { // for each fleet...
      arma::vec effort_now = Rcpp::as<arma::vec>(fleet_fishing_effort[FLEET]); 
      master_effort_by_fleet.slice(FLEET).col(MONTH) = effort_now; // format (cell, month, fleet) 
    } // record how much fishing there will be this month and where.
    
    Rcpp::Rcout << "fishers at the ready - ";

    
    next_pop.zeros(); // start with a fresh object, which we will fill for this month, to give to the next month.
    for (int FLEET = 0; FLEET < n_fleets; FLEET++) {catch_weight[FLEET] = arma::vec(max_cell, arma::fill::zeros);} // reset the catch_weight this month, otherwise you're accumulating catch across 12 months

    
    for (int AGE = 0; AGE < max_age; AGE++) {
      
      arma::mat current_pop_AGE = current_pop.slice(AGE); // obtain the current age's population once
      
      
    // 3. move fish -----------------------------------------------------------

      arma::mat moved_population = movement_function(AGE, // for every aged fish, move them to another spot
                                                     max_cell,
                                                     adult_movement_prob,
                                                     current_pop_AGE);

      current_pop.slice(AGE) = moved_population; // cell x length x age

    
    // 4. kill fish -----------------------------------------------------------

      Rcpp::List mortality_outputs = mortality_function(AGE, // for every aged fish, they are exposed to natural and fishing mortalities
                                                        MONTH,
                                                        max_cell,
                                                        n_lengths,
                                                        natural_mortality,
                                                        weight, // vec of size n_lengths
                                                        sel_this_year, // vec n_lengths
                                                        current_pop_AGE, // cell x length
                                                        fleet_fishing_effort, // from the effort_function - list of size n_fleets, made of vectors of size n_cell,
                                                        catchability_by_fleet
      );
      arma::mat survived = Rcpp::as<arma::mat>(mortality_outputs["tot_survived"]);
      
      // also fill the master outputs
      master_pop_survived(MONTH).slice(AGE) = survived;
      
      arma::mat F_this_month = Rcpp::as<arma::mat>(mortality_outputs["fishing_mortality"]); // cell x length

      master_F_numerator   += F_this_month % current_pop_AGE;  // weight F by abundance
      master_F_denominator += current_pop_AGE;                 // accumulate weights
      
      
      // extract the catch
      Rcpp::List catch_weights_this_age = Rcpp::as<Rcpp::List>(mortality_outputs["catch_weight_by_fleet"]);
      Rcpp::List catch_numbers_this_age = Rcpp::as<Rcpp::List>(mortality_outputs["catch_numbers_by_fleet"]);
      
      for (int FLEET = 0; FLEET < n_fleets; FLEET++) { // for each fleet...
        
        arma::vec accum = Rcpp::as<arma::vec>(catch_weight[FLEET]);
        accum += Rcpp::as<arma::vec>(catch_weights_this_age[FLEET]);
        catch_weight[FLEET] = accum; // update the catch so that next month, the effort can be distributed accordingly
        
        // and record it
        master_catch_weight_total(FLEET).col(AGE) = Rcpp::as<arma::vec>(catch_weights_this_age[FLEET]);
        master_catch_number_total(FLEET).slice(AGE) = Rcpp::as<arma::mat>(catch_numbers_this_age[FLEET]);
      } // record how much catch there is this month and where.
      
      
      // 5. grow fish ---------------------------------------------------------
      
      if (MONTH < 11) { // Jan-Nov, fish grow in length, but don't age.
          next_pop.slice(AGE) = survived * age_transition; // max_cell x n_lengths  
        } 
      else if (MONTH == 11) { // in December, fish grow AND age.
        if (AGE < (max_age - 1)) {
          next_pop.slice(AGE + 1) = survived * age_transition; // max_cell x n_lengths
        }
        // the oldest fish (AGE = max_age) are discarded and not tracked in the model anymore
      } // note: the age-transition matrix MUST be calculated on a monthly basis.
      
    } // end of AGE loop
    
    Rcpp::Rcout << "fish successfully moved, killed, and grown - ";
    
    
    // 6. recruit fish --------------------------------------------------------

    if (is_spawn_month[MONTH]) { // if it's spawning month...

      Rcpp::List recruitment_outputs = recruitment_function(MONTH,          // current spawning month
                                                            max_cell,       // number of cells in the grid
                                                            max_age,        // oldest age of the species
                                                            BHa,            // Beverton-Holt parameter 
                                                            BHb,            // Beverton-Holt parameter
                                                            PF,             // proportion of females 
                                                            ha_scaling,     // hyperallometry scaling factor
                                                            maturity,       // vec of size ncell of the proportion of the each length to be mature.
                                                            weight,         // vec of size ncell of the weight of each length.
                                                            settlement,     // vector of size max_cell giving probability of recruiting in each cell (based on habitat etc.)
                                                            current_pop     // numbers of fish in each cell x length x age
      );

      arma::vec settle_recs = recruitment_outputs["settle_recs"];
      january_recruits += settle_recs; // sum the recruits of all spawning months
      Rcpp::Rcout << "Spawning MONTH=" << MONTH + 1
                  << " sum(january_recruits)=" << arma::sum(january_recruits)
                  << std::endl;
      
    } // end of spawning MONTH loop
    
    Rcpp::Rcout << "fish spawned." << std::endl;
    
    
    // 7. handoff the processed population to next month ----------------------
    
    if (MONTH < 11) {
      current_pop = next_pop;
      
      // age structure: sum all cells and lengths for each age slice
      Rcpp::Rcout << "age structure: ";
      for (int a = 0; a < next_pop.n_slices; a++) {
        Rcpp::Rcout << arma::accu(next_pop.slice(a)) << " ";
      }
      Rcpp::Rcout << std::endl;
      
      // size structure: sum all cells and ages for each length
      Rcpp::Rcout << "size structure: ";
      arma::mat summed_ages = arma::zeros(next_pop.n_rows, next_pop.n_cols);
      for (int a = 0; a < next_pop.n_slices; a++) {
        summed_ages += next_pop.slice(a);
      }
      arma::rowvec size_structure = arma::sum(summed_ages, 0); // sum over cells (rows)
      Rcpp::Rcout << size_structure << std::endl;
    }
    
  } // end of the MONTH loop

  
  // compute weighted mean annual F
  master_F_denominator.replace(0, 1e-10);
  master_F = master_F_numerator / master_F_denominator; // cell x length
 
  // compute SSB per spawning month
  arma::vec ha_weight = arma::pow(weight, ha_scaling);
  arma::vec maturity_x_haweight = maturity % ha_weight;
  arma::vec master_SSB_total(max_cell, arma::fill::zeros); // summed across spawning months
 
  for (arma::uword k = 0; k < spawning_months.n_elem; k++) {
    int sp_month = static_cast<int>(spawning_months(k)) - 1; // zero-indexed
    arma::cube pop_this_month = master_pop_survived(sp_month); // cell x length x age
    for (int AGE = 0; AGE < max_age; AGE++) {
      master_SSB_total += PF * (pop_this_month.slice(AGE) * maturity_x_haweight);
    }
  }
 
  next_pop.slice(0).col(0) = january_recruits; // ← your existing line
  
  Rcpp::Rcout << "fish successfully recruited, population handed off to next month." << std::endl;
  
  
  return Rcpp::List::create(
    
    // also add master outputs
    Rcpp::Named("current_pop") = current_pop, // outputs - 
    Rcpp::Named("next_pop") = next_pop, // outputs - 
    
    Rcpp::Named("annual_F")              = master_F,          // cell x length
    Rcpp::Named("spawning_biomass")      = master_SSB_total,  // cell, summed across spawning months
    
    Rcpp::Named("catch_number_by_fleet") = master_catch_number_total, // cell x month x fleet
    Rcpp::Named("catch_weight_by_fleet") = master_catch_weight_total, // cell x month x fleet
    Rcpp::Named("master_pop_survived") = master_pop_survived // cell x month x age
  
  );
}
