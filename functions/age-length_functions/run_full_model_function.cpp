#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include "internals/effort_function.h"
#include "internals/mortality_function.h"
#include "internals/pseudo_effort_function_2.h"
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
    arma::mat juv_movement_prob, // max_cell x max_cell
    arma::mat spawn_movement_prob, // max_cell x max_cell
    
    Rcpp::CharacterVector fleet_names,
    Rcpp::List fleet_info
) {
  
  
  // 0a. create master objects ------------------------------------------------
  
  // 1 where is the fish every month x year 
  arma::cube next_pop(max_cell, n_lengths, max_age, arma::fill::zeros); // function object, to put survivors and recruits in 

  arma::field<arma::cube> master_current_pop(12); // format month x (cell x length x age)
  for (int m = 0; m < 12; m++) {
    master_current_pop(m) = arma::cube(max_cell, n_lengths, max_age, arma::fill::zeros);
  }  
  
  // 2 where is the effort every month x year
  int n_fleets = fleet_names.size();
  Rcpp::List fleet_fishing_effort(n_fleets); // function object
  arma::cube master_effort_by_fleet(max_cell, 12, n_fleets, arma::fill::zeros); // format cell x month x fleet
  
  // 3 how much catch is there every year
  arma::vec master_yearly_catch(n_fleets, arma::fill::zeros); // vec of length n_fleet
  arma::cube master_catch_by_length(max_cell, n_lengths, n_fleets, arma::fill::zeros); // cell x lengths x fleet, this is an important output to understand whether shore fishers are catching less fish over time
    
  // 4 spawning biomass every year
  arma::vec january_recruits(max_cell, arma::fill::zeros); // function object, sum of recruits of all spawning months. gets renewed every year.
  double master_SSB = 0.0; // single value (should be the mean SSB across spawning months)
  
  // 5 fishing mortality
  arma::mat master_F(12, n_fleets, arma::fill::zeros); 
  
  
  // 0b. extract objects that need to be calculated only once per year --------
  
  arma::vec sel_this_year = selectivity.col(YEAR);
  
  std::vector<arma::mat> catchability_by_fleet(n_fleets);
  for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
    Rcpp::List fl = fleet_info[FLEET];
    catchability_by_fleet[FLEET] = Rcpp::as<arma::cube>(fl["catchability"]).slice(YEAR); // becomes a list of fleet x (max_cell x 12)
    }
  
  std::vector<bool> is_spawn_month(12, false);
  for (arma::uword k = 0; k < spawning_months.n_elem; k++)
    is_spawn_month[static_cast<int>(spawning_months(k)) - 1] = true;
  arma::vec ha_weight = arma::pow(weight, ha_scaling);
  arma::vec maturity_x_haweight = maturity % ha_weight;
  
  
  // START OF THE FUNCTION HERE ===============================================
  
  for (int MONTH = 0; MONTH < 12; MONTH++) { 
    
    Rcpp::Rcout << "Year: " << YEAR+1 << ", Month: " << MONTH+1 << std::endl;
    
    // 1. figure out where high CPUE was in the previous month ------------------
    
    Rcpp::List expected_catch(n_fleets); // (list n_fleet, obj vec ncell)

    if (YEAR == 0 && MONTH == 0) { // for the first month of the first year, no information about catch, so everywhere is the same.
      for (int FLEET = 0; FLEET < fleet_names.size(); FLEET++) {
        expected_catch[FLEET] = arma::vec(max_cell, arma::fill::ones); // (list n_fleet, obj vec ncell)

      } 
    } // for all fleets for year1 month1, fill the expected catch (and expected catch squared) with zeroes, since we don't have any previous month CPUE.
    
    else {
      
      Rcpp::List pseudo = pseudo_effort_function(MONTH,
                                                 max_cell,
                                                 n_lengths,
                                                 max_age,
                                                 natural_mortality,
                                                 sel_this_year, // selectivity-retention of fish in this year (extracted in the master function). currently a single selectivity for all fleets
                                                 current_pop, // cell x length
                                                 catchability_by_fleet // extracted before the MONTH loop in the master function
      );
      
      expected_catch = pseudo["pseudo_catch"]; // (list n_fleet, obj vec ncell)

    } // for all months and years after y1 m1, use the previous month's catch and effort distribution to obtain pseudo catch for this month
    
    Rcpp::Rcout << "fisher strategy done - ";
    
    
    // 2. distribute effort accordingly ---------------------------------------
    
    Rcpp::List effort_distribution_output = distribute_effort_function(MONTH,
                                                                       YEAR,
                                                                       max_cell,
                                                                       expected_catch,
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


      // 3. move fish -----------------------------------------------------------
      
      for (int AGE = 0; AGE < max_age; AGE++) {

      arma::mat moved_population = movement_function(AGE, // for every aged fish, move them to another spot
                                                     MONTH, 
                                                     max_cell,
                                                     adult_movement_prob, // ncell x ncell
                                                     juv_movement_prob, // ncell x ncell
                                                     spawn_movement_prob, // ncell x ncell
                                                     current_pop.slice(AGE)); // ncell x nlengths

      current_pop.slice(AGE) = moved_population; // cell x length x age

      } // end of AGE loop for movement

      Rcpp::Rcout << "fish successfully moved - ";

      
      // 4. kill fish -----------------------------------------------------------

      for (int AGE = 0; AGE < max_age; AGE++) {
      
      Rcpp::List mortality_outputs = mortality_function(AGE, // for every aged fish, they are exposed to natural and fishing mortalities
                                                        MONTH,
                                                        max_cell,
                                                        n_lengths,
                                                        natural_mortality,
                                                        weight, // vec of size n_lengths
                                                        sel_this_year, // vec n_lengths
                                                        current_pop.slice(AGE), // cell x length
                                                        fleet_fishing_effort, // from the effort_function - list of size n_fleets, made of vectors of size n_cell,
                                                        catchability_by_fleet);
      arma::mat survived = Rcpp::as<arma::mat>(mortality_outputs["tot_survived"]);

      // extract the catch
      Rcpp::List catch_weights_this_age = Rcpp::as<Rcpp::List>(mortality_outputs["catch_weight_by_fleet"]); // list of n_fleets x (cell x lengths)
      Rcpp::List catch_numbers_this_age = Rcpp::as<Rcpp::List>(mortality_outputs["catch_numbers_by_fleet"]); // list of n_fleets x (cell x lengths)
      
      // and the fishing mortality
      arma::cube F_by_fleet_this_age = Rcpp::as<arma::cube>(mortality_outputs["fishing_mortality_by_fleet"]); // cell x length x fleet
      for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
        master_F(MONTH, FLEET) += arma::accu(F_by_fleet_this_age.slice(FLEET));
      }
      
      for (int FLEET = 0; FLEET < n_fleets; FLEET++) { // for each fleet...
        
        master_yearly_catch[FLEET] += arma::accu(Rcpp::as<arma::vec>(catch_weights_this_age[FLEET])); // accumulate the catch of every month, and sum over cells.
        
        master_catch_by_length.slice(FLEET) += Rcpp::as<arma::mat>(catch_numbers_this_age[FLEET]); // accumulate the catch of every month, but keep the spatial structure
        
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
      
    } // end of AGE loop for mortality + growth
    
    // --- CHECK catch of each fleet ---

    Rcpp::Rcout << " total pop: " << arma::accu(next_pop) << std::endl;
    Rcpp::Rcout << " | total catch: " ;
    for (int FLEET = 0; FLEET < n_fleets; FLEET++) {
      Rcpp::Rcout << master_yearly_catch[FLEET] << " ";
      }
    
    Rcpp::Rcout << "fish successfully killed, and grown - " << std::endl;
    
    
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
      
      double ssb_this_month = 0.0;
      for (int AGE = 0; AGE < max_age; AGE++) {
        ssb_this_month += arma::accu(next_pop.slice(AGE) * maturity_x_haweight);
      }
      master_SSB += PF * ssb_this_month / spawning_months.n_elem;      
    
      Rcpp::Rcout << "fish spawned." << std::endl;
    
    } // end of spawning MONTH loop
    
    
    // 7. handoff the processed population to next month ----------------------
    
    current_pop = next_pop;
    master_current_pop(MONTH) = current_pop; // cell x length x age snapshot for this month
    
    
  } // end of the MONTH loop
  
  
  next_pop.slice(0).col(0) = january_recruits; // assign recruits to the youngest age & smallest length of the next year
  
  Rcpp::Rcout << "fish successfully recruited, population handed off to next month." << std::endl;
  
  
  return Rcpp::List::create(
    Rcpp::Named("next_pop")           = next_pop, // cube (cell x lengt), needed to hand the december population to January of the next year in R 
    Rcpp::Named("master_current_pop") = master_current_pop, // month x (cell x length x age)

    Rcpp::Named("effort_by_fleet")    = master_effort_by_fleet, // cell x month x fleet
    
    Rcpp::Named("yearly_catch")       = master_yearly_catch, // vec of size n_fleet
    Rcpp::Named("catch_by_length")    = master_catch_by_length,
    
    Rcpp::Named("fishing_mortality")  = master_F, // month x fleet
    
    Rcpp::Named("spawning_biomass")   = master_SSB  // single value
  );
}
