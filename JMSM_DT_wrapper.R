# Wrapper Function to implement JMSM DT method in paper
JMSM_DT <- function(data, treatment, time_varying_cov, interval){
  
  print("weight estimaion for A1")
  weight_JMSM_DT_A1 <- weight_JMSM_DT(data = simulation_data_ragged,
                                      treatment = "A1",interval = interval)
  print("weight estimaion for A2")
  weight_JMSM_DT_A2 <- weight_JMSM_DT(data = simulation_data_ragged,
                                      treatment = "A2",interval = interval)
  print("marginal structural model")
  marginal_structural_cox_bias_result <-
    marginal_structural_cox_bias(data = simulation_data_ragged,
                                 weight = matrix(
                                   c(weight_JMSM_DT_A1, weight_JMSM_DT_A2),
                                   nrow = 2,
                                   byrow = T
                                 ))
  return(marginal_structural_cox_bias_result)
}
