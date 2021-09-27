# Function to fit the marginal structural cox model using the estimated weight and calculate the bias for the phis
marginal_structural_cox_bias <- function(data, weight){
  # First aggregate the data into patient level
  data_aggregated <- data %>% 
    group_by(id) %>% 
    summarise(A1 = sum(A1),
              A2 = sum(A2)) %>% 
    mutate(A1 = ifelse(A1 >= 1, 1, A1),
           A2 = ifelse(A2 >= 1, 1, A2)) %>% 
    inner_join(data %>% 
                 filter(Y == 1) %>% 
                 group_by(id) %>%
                 filter(row_number()== 1) %>% 
                 ungroup %>% 
                 bind_rows(tibble(id = setdiff(1:1000, data %>% 
                                                 filter(Y == 1) %>% pull(id)),
                                  time = 49.5,
                                  Y = 0)) %>% 
                 select(id, time, Y))
  
  marginal_structural_survival_model_cox_weight <- data_aggregated %>% 
    mutate(weight = apply(weight,2, prod)) %>% 
    coxph(Surv(time, Y) ~ A1 + A2, data = .)
  cox_phi_bias <- marginal_structural_survival_model_cox_weight %>% 
    broom::tidy() %>% 
    mutate(bias = (estimate - c(-0.5, -0.3))/c(6,-60)) %>% 
    select(term, bias)
  return(cox_phi_bias)
}