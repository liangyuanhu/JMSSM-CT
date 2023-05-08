##############################################################################################################################################################################################
# R codes for the simulation part of the paper:
# Estimating the causal effects of multiple longitudinal treatments in continuous time with applications to COVID-19
#####################################################################################################################################################################################################################

# Source functions to estimate weight using time-varying Cox model
source("functions/weight_Cox_TV_wrapper.R")
source("functions/weight_Cox_TV.R")

# Source functins to estimate weight using time-varying RRF
source("functions/weight_RRF_TV_forest_wrapper.R")
source("functions/weight_RRF_TV_forest.R")

# Source functions to estimate weight using discrete time random forest
source("functions/weight_DT_RF_wrapper.R")
source("functions/weight_DT_RF.R")

# Source wrapper function to implement JMSM DT
source("functions/JMSM_DT_wrapper.R")

# Source wrapper function to implement JMSSM CT
source("functions/JMSSM_CT_wrapper.R")

# Source function to implement marginal structural cox
source("functions/marginal_structural_cox.R")

# Source data generation function
source("functions/sim_marginal_structural_cox.R")

# Source other supportive functions 
source("functions/reformat_to_recurrent.R")
source("functions/zero_inflated_Poisson.R")


library(LTRCforests)
library(survival)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(survival) 
library(randomForest)
##############################################################################
# Simulate the data #####################################################################################################################################################################################################################
set.seed(11112)
simulation_data_aligned <- sim_marginal_structural_cox(subjects = 1000, tpoints = 100, psiA1 = -.5, psiA2 = -.3, n = 1)
# simulation_data_aligned <- sim_marginal_structural_cox(subjects = 500, tpoints = 100, psiA1 = -.5, psiA2 = -.3, n = 1)
# simulation_data_aligned <- sim_marginal_structural_cox(subjects = 250, tpoints = 100, psiA1 = -.5, psiA2 = -.3, n = 1)
# Randomly select a varying number of observations for each individual to create a “ragged” longitudinal dataset
simulation_data_ragged <- simulation_data_aligned %>% 
  filter(time != 0) %>% 
  sample_frac(0.7) %>% 
  ungroup() %>% 
  bind_rows(simulation_data_aligned %>% 
              filter(time == 0)) %>% 
  arrange(id,time) 

##############################################################################
# Test the 2 main functions for the 3 methods:  ################################################
# 1) JMSSM CT Cox, 2) JMSSM CT RRF-TV forest, 3) JMSM DT  ################################################
#####################################################################################################################################################################################################################

# JMSSM CT Cox: Weight estimated by Cox, smoothing the Nelson-Aalen estimator----
set.seed(111)
result_JMSSM_CT_Cox_weight <- JMSSM_CT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), smooth = T, weight_method = "Cox")
# [1] "weight estimaion for A1"
# [1] "weight estimaion for A2"
# [1] "marginal structural model"
result_JMSSM_CT_Cox_weight
# Test the function
# term    bias
# <chr>  <dbl>
# 1 A1    0.0833
# 2 A2    0.0947 

# JMSSM CT RRF-TV forest: Weight estimated by RRF-TV forest, smoothing the Nelson-Aalen estimator -----------------
set.seed(111)
result_JMSSM_CT_RRF_TV_forest_weight <- JMSSM_CT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), smooth = T, weight_method = "RRF-TV forest")
# [1] "weight estimaion for A1"
# [1] "weight estimaion for A2"
# [1] "marginal structural model"
# Test the function 
result_JMSSM_CT_RRF_TV_forest_weight
# term  bias
# <chr>    <dbl>     
# 1 A1      0.032   
# 2 A2      0.017   

# JMSM DT: Weight estimated by discrete time random forest--------
result_JMSM_DT <- JMSM_DT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), interval = 2)
result_JMSM_DT
# Test the function
# term    bias
# <chr>  <dbl>
# 1 A1    0.0333
# 2 A2    0.0307 
##############################################################################
##############################################################################
# Implement simulation: part 1  #####################################################################################################################################################################################################################
#####################################
# 250 repetitions for method (i), (ii), (iii), (iv)
result_250_repetitions_JMSSM_CT_Cox_weight_no_smooth <- NULL
result_250_repetitions_JMSSM_CT_Cox_weight_with_smooth <- NULL
result_250_repetitions_JMSSM_CT_RRF_TV_forest_weight_no_smooth <- NULL
result_250_repetitions_JMSSM_CT_RRF_TV_forest_weight_with_smooth <- NULL

for (i in 1:250){
  set.seed(11112)
  simulation_data_aligned <- sim_marginal_structural_cox(subjects = 1000, tpoints = 100, psiA1 = -.5, psiA2 = -.3, n = 1)
  # simulation_data_aligned <- sim_marginal_structural_cox(subjects = 500, tpoints = 100, psiA1 = -.5, psiA2 = -.3, n = 1)
  # simulation_data_aligned <- sim_marginal_structural_cox(subjects = 250, tpoints = 100, psiA1 = -.5, psiA2 = -.3, n = 1)
  # Randomly select a varying number of observations for each individual to create a “ragged” longitudinal dataset
  simulation_data_ragged <- simulation_data_aligned %>% 
    filter(time != 0) %>% 
    sample_frac(0.7) %>% 
    ungroup() %>% 
    bind_rows(simulation_data_aligned %>% 
                filter(time == 0)) %>% 
    arrange(id,time) 
  
  # (i) JMSSM CT Cox: Weight estimated by Cox, without smoothing the Nelson-Aalen estimator----
  result_JMSSM_CT_Cox_weight_no_smooth <- JMSSM_CT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), smooth = F, weight_method = "Cox")
  
  # (ii) JMSSM CT Cox: Weight estimated by Cox, smoothing the Nelson-Aalen estimator----
  result_JMSSM_CT_Cox_weight_with_smooth <- JMSSM_CT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), smooth = T, weight_method = "Cox")
  
  # (iii) JMSSM CT RRF-TV forest: Weight estimated by RRF-TV forest, without smoothing the Nelson-Aalen estimator -----------------
  result_JMSSM_CT_RRF_TV_forest_weight_no_smooth <- JMSSM_CT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), smooth = F, weight_method = "RRF-TV forest")
  
  # (iv) JMSSM CT RRF-TV forest: Weight estimated by RRF-TV forest, without smoothing the Nelson-Aalen estimator -----------------
  result_JMSSM_CT_RRF_TV_forest_weight_with_smooth <- JMSSM_CT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), smooth = T, weight_method = "RRF-TV forest")
  result_250_repetitions_JMSSM_CT_Cox_weight_no_smooth <- result_250_repetitions_JMSSM_CT_Cox_weight_no_smooth %>% bind_cols(result_JMSSM_CT_Cox_weight_no_smooth)
  result_250_repetitions_JMSSM_CT_Cox_weight_with_smooth <- result_250_repetitions_JMSSM_CT_Cox_weight_no_smooth %>% bind_cols(result_JMSSM_CT_Cox_weight_with_smooth)
  result_250_repetitions_JMSSM_CT_RRF_TV_forest_weight_no_smooth <- result_250_repetitions_JMSSM_CT_RRF_TV_forest_weight_no_smooth %>% bind_cols(result_JMSSM_CT_RRF_TV_forest_weight_no_smooth)
  result_250_repetitions_JMSSM_CT_RRF_TV_forest_weight_with_smooth <- result_250_repetitions_JMSSM_CT_RRF_TV_forest_weight_with_smooth %>% bind_cols(result_JMSSM_CT_RRF_TV_forest_weight_with_smooth)
  print("Finish", i)
}
##############################################################################
##############################################################################
# Implement simulation: part 2  #####################################################################################################################################################################################################################

# 250 repetitions for comparing method:JMSM-DT vs JMSSM-CT
result_250_repetitions_JMSSM_CT_RRF_TV_forest_weight_with_smooth <- NULL
result_250_repetitions_JMSM_DT_weight_2days_interval <- NULL
result_250_repetitions_JMSM_DT_weight_1days_interval <- NULL
result_250_repetitions_JMSM_DT_weight_0.5days_interval <- NULL

for (i in 1:250){
  result_JMSSM_CT_RRF_TV_forest_weight_with_smooth <- JMSSM_CT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), smooth = T, weight_method = "RRF-TV forest")
  result_JMSM_DT_weight_2days_interval <- MSM_DT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), interval = 2)
  result_JMSM_DT_weight_1days_interval <- MSM_DT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), interval = 1)
  result_JMSM_DT_weight_0.5days_interval <- MSM_DT(data = simulation_data_ragged, treatment= c("A1", "A2"), time_varying_cov = c("L1", "L2"), interval = 0.5)
  result_250_repetitions_JMSSM_CT_RRF_TV_forest_weight_with_smooth <- result_250_repetitions_JMSSM_CT_RRF_TV_forest_weight_with_smooth %>% bind_cols(result_JMSSM_CT_RRF_TV_forest_weight_with_smooth)
  result_250_repetitions_JMSM_DT_weight_2days_interval <- result_250_repetitions_JMSM_DT_weight_2days_interval %>% bind_cols(result_JMSM_DT_weight_2days_interval)
  result_250_repetitions_JMSM_DT_weight_1days_interval <- result_250_repetitions_JMSM_DT_weight_2days_interval %>% bind_cols(result_JMSM_DT_weight_1days_interval)
  result_250_repetitions_JMSM_DT_weight_0.5days_interval <- result_250_repetitions_JMSM_DT_weight_2days_interval %>% bind_cols(result_JMSM_DT_weight_0.5days_interval)
  print("Finish", i)
}
