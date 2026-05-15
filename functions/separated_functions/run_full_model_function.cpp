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
    const int max_year,
    
    arma::cube yearly_pop,
    arma::mat weight,
    arma::cube selectivity,
    const double natural_mortality,
    
    arma::vec spawning_months,
    const double BHa,
    const double BHb,
    const double PF,
    const double ha_scaling,
    arma::mat maturity, 
    arma::vec settlement,
    
    arma::mat adult_movement_prob,
    
    Rcpp::CharacterVector fleet_names,
    Rcpp::List fleet_info
) {
  
  
  // 0. create objects to store the master outputs ----------------------------
  
  int n_fleets = fleet_names.size();
  
  arma::cube              master_pop_survived(max_cell, 12, max_age, arma::fill::zeros); // format (cell, month, age) 
  arma::cube              master_effort_by_fleet(max_cell, 12, n_fleets, arma::fill::zeros); // format (cell, month, fleet) 
  arma::field<arma::cube> master_catch_number_total(n_fleets); // format cell x month x age x fleet
  arma::field<arma::cube> master_catch_weight_total(n_fleets); // format cell x month x age x fleet
  
  arma::vec               january_recruits(max_cell, arma::fill::zeros); // this is for a sum of recruits of all spawning months. gets renewed every year.
  Rcpp::List              catch_weight(n_fleets);
  Rcpp::List              fleet_fishing_effort(n_fleets);
  
  arma::mat               F(max_cell, max_age, arma::fill::zeros); // just to test whether fishing mortality is realistic
  
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    master_catch_number_total(FLEET) = arma::cube(max_cell, 12, max_age, arma::fill::zeros);
    master_catch_weight_total(FLEET) = arma::cube(max_cell, 12, max_age, arma::fill::zeros);
    
    catch_weight[FLEET] = arma::vec(max_cell, arma::fill::ones);
    fleet_fishing_effort[FLEET]  = arma::vec(max_cell, arma::fill::ones);
  } // fill the master outputs with zeroes to start with, so that they're the right dimensions
  
  
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
    
    // 2. move fish -------------------------------------------------------------
    
    for (int AGE = 0; AGE < max_age; AGE++) { 
      
      arma::vec moved_population = movement_function(AGE,
                                                     MONTH, 
                                                     max_cell, 
                                                     adult_movement_prob, 
                                                     yearly_pop);
      
      yearly_pop.slice(AGE).col(MONTH) = moved_population; // cell x month x age

      } // for every aged fish, move them to another spot
    
    
    // 3. distribute effort -----------------------------------------------------

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
      master_effort_by_fleet.slice(FLEET).col(MONTH) = Rcpp::as<arma::vec>(fleet_fishing_effort[FLEET]); // format (cell, month, fleet) 
    } // record how much fishing there is this month and where.
    
    
    
    // 4. kill fish -------------------------------------------------------------
    
    for (int FLEET = 0; FLEET < n_fleets; FLEET++) { // reset the catch_weight
      catch_weight[FLEET] = arma::vec(max_cell, arma::fill::zeros);
    }

    for (int AGE = 0; AGE < max_age; AGE++) {
      
      Rcpp::List mortality_outputs = mortality_function(AGE,
                                                        MONTH,
                                                        YEAR,
                                                        max_cell,
                                                        natural_mortality,
                                                        weight, // age x month
                                                        selectivity, // age x month x year - likely will need to change by fleet
                                                        yearly_pop, // cell x month x age
                                                        fleet_fishing_effort, // from the effort_function - list of size n_fleets, made of vectors of size n_cell,
                                                        fleet_info
      );
      
      // also fill the master outputs
      arma::vec survived = Rcpp::as<arma::vec>(mortality_outputs["tot_survived"]);
      
      master_pop_survived.slice(AGE).col(MONTH) = survived;
      
      F.col(AGE) = Rcpp::as<arma::vec>(mortality_outputs["fishing_mortality"]);
      
      // extract the catch
      Rcpp::List catch_weights_this_age = Rcpp::as<Rcpp::List>(mortality_outputs["catch_weight_by_fleet"]);
      Rcpp::List catch_numbers_this_age = Rcpp::as<Rcpp::List>(mortality_outputs["catch_numbers_by_fleet"]);
      
      for (int FLEET = 0; FLEET < n_fleets; FLEET++) { // for each fleet...
        
        arma::vec accum = Rcpp::as<arma::vec>(catch_weight[FLEET]);
        accum += Rcpp::as<arma::vec>(catch_weights_this_age[FLEET]);
        catch_weight[FLEET] = accum; // update the catch so that next month, the effort can be distributed accordingly
        
        // and record it
        master_catch_weight_total(FLEET).slice(AGE).col(MONTH) = Rcpp::as<arma::vec>(catch_weights_this_age[FLEET]);
        master_catch_number_total(FLEET).slice(AGE).col(MONTH) = Rcpp::as<arma::vec>(catch_numbers_this_age[FLEET]);
      } // record how much catch there is this month and where.
      
      
      if (MONTH < 11) {
        if (AGE < (max_age - 1)) {
          yearly_pop.slice(AGE).col(MONTH + 1) = survived;
        }
      } else if (MONTH == 11) {
        if (AGE < (max_age - 1)) {
          yearly_pop.slice(AGE + 1).col(0) = survived;
        } 
      } // in december, age-up the fish, or assign this population to the next month for all other months. skip the oldest fish, who aren't taken into account anymore.
      
    }

    
    // 5. recruit fish ----------------------------------------------------------

    arma::uvec spawn_check = arma::find(spawning_months - 1 == MONTH);

    if (spawn_check.n_elem > 0) { // if it's spawning month...

      Rcpp::List recruitment_outputs = recruitment_function(MONTH,          // current spawning month
                                                            max_cell,       // number of cells in the grid
                                                            max_age,        // oldest age of the species
                                                            BHa,            // Beverton-Holt parameter 
                                                            BHb,            // Beverton-Holt parameter
                                                            PF,             // proportion of females 
                                                            ha_scaling,     // hyperallometry scaling factor
                                                            maturity,       // age x month matrix of the proportion of the population to be mature over 12 months.
                                                            weight,         // age x month matrix of the weight of a fish of each age over months over 12 months.
                                                            settlement,     // vector of size max_cell giving probability of recruiting in each cell (based on habitat etc.)
                                                            yearly_pop      // numbers of fish in each cell x month x age
      );

      arma::vec settle_recs = recruitment_outputs["settle_recs"];
      january_recruits += settle_recs; // sum the recruits of all spawning months
      Rcpp::Rcout << "Spawning MONTH=" << MONTH 
                  << " sum(settle_recs)=" << arma::sum(settle_recs)
                  << " sum(january_recruits)=" << arma::sum(january_recruits)
                  << std::endl;
    }
    
  } // end of the MONTH loop
  
  yearly_pop.slice(0).col(0) = january_recruits; // add recruits to january
  
  return Rcpp::List::create(
    
    // also add master outputs
    Rcpp::Named("yearly_pop") = yearly_pop, // outputs - 
    Rcpp::Named("fishing_mortalities") = F, // checks only - 
    
    Rcpp::Named("catch_number_by_fleet") = master_catch_number_total, // cell x month x fleet
    Rcpp::Named("catch_weight_by_fleet") = master_catch_weight_total, // cell x month x fleet
    Rcpp::Named("master_pop_survived") = master_pop_survived // cell x month x age
  
  );
}
