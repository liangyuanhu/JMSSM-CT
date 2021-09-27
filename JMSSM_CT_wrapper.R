# Wrapper Function to implement JMSM CT method in paper

JMSSM_CT <- function(data, treatment, time_varying_cov, smooth, weight_method){
  if (weight_method == "Cox"){
    print("weight estimaion for A1")
    weight_Cox_A1 <- weight_Cox_treatment(data = simulation_data_ragged, treatment = "A1", time_varying_cov = c("L1", "L2"), smooth = smooth)
    print("weight estimaion for A2")
    weight_Cox_A2 <- weight_Cox_treatment(data = simulation_data_ragged, treatment = "A2", time_varying_cov = c("L1", "L2"), smooth = smooth)
    print("marginal structural model")
    marginal_structural_cox_bias_result <-
      marginal_structural_cox_bias(data = simulation_data_ragged,
                                   weight = matrix(
                                     c(weight_Cox_A1, weight_Cox_A2),
                                     nrow = 2,
                                     byrow = T
                                   ))
    
  }
  if (weight_method == "RRF-TV forest"){
    print("weight estimaion for A1")
    weight_LTRCforests_treatment_A1 <- weight_LTRCforests_treatment(data = simulation_data_ragged, treatment = "A1", time_varying_cov = c("L1", "L2"))
    print("weight estimaion for A2")
    weight_LTRCforests_treatment_A2 <- weight_LTRCforests_treatment(data = simulation_data_ragged, treatment = "A2", time_varying_cov = c("L1", "L2"))
    print("marginal structural model")
    marginal_structural_cox_bias_result <-
      marginal_structural_cox_bias(data = simulation_data_ragged,
                                   weight = matrix(
                                     c(weight_LTRCforests_treatment_A1, weight_LTRCforests_treatment_A2),
                                     nrow = 2,
                                     byrow = T
                                   ))
  }
  return(marginal_structural_cox_bias_result)
}
