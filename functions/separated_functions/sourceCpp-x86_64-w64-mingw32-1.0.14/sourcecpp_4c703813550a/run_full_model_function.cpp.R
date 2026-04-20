`.sourceCpp_1_DLLInfo` <- dyn.load('C:/Users/localuser/Desktop/GitHub/Fournier-Carnoy_Yijarup_Wadandi_population_model/functions/separated_functions/sourceCpp-x86_64-w64-mingw32-1.0.14/sourcecpp_4c703813550a/sourceCpp_4.dll')

test_wrapper <- Rcpp:::sourceCppFunction(function(x) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_test_wrapper')
run_full_model_function <- Rcpp:::sourceCppFunction(function(YEAR, max_cell, max_age, max_year, yearly_pop, weight, selectivity, natural_mortality, spawning_months, BHa, BHb, PF, ha_scaling, maturity, settlement, adult_movement_prob, fleet_names, fleet_info) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_run_full_model_function')

rm(`.sourceCpp_1_DLLInfo`)
