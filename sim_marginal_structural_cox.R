sim_marginal_structural_cox <- function(subjects = 2500, tpoints = 10, psiA1 = .5, psiA2 = .3, n = 1){
  cat("Datasets being saved in:", getwd(), "\n")
  for (select.seed.val  in 1:n){
    # seeds
    #beta =parameters for S and L
    beta0 = log(3/7)
    beta1 = 2
    beta2 = log(0.5)
    beta3 = log(1.5)
    #alphaA1_ = parameters for A1
    alphaA1_0= -log(2/7)
    alphaA1_1= -0.5
    alphaA1_2 = -0.5
    alphaA1_3 = -log(3/5)
    alphaA1_4 = -1.2
    alphaA1_22 = 0.5
    alphaA1_32 = 1.1
    alphaA1_42 = -0.8
    alphaA2_0= -log(2/7)
    alphaA2_1= -0.5
    alphaA2_2 = -0.5
    alphaA2_3 = -log(3/5)
    alphaA2_4 = -1.2
    alphaA2_22 = 0.5
    alphaA2_32 = 1.1
    alphaA2_42 = -0.8
    cut=30;
    lambda0=0.005 #scale parameter for weibull
    shape=1 #shape parameter for weibull
    
    outputx <- NULL
    set.seed(select.seed.val)
    
    pb <- txtProgressBar(min = 0, max = subjects, style = 3)
    
    for (id in 1:subjects){
      expvar = rexp(1, rate = 1)
      T0 = (expvar/(lambda0^shape))^(1/shape)
      if (T0<cut) {IT0=1} else {IT0=0}
      maxT = tpoints
      initA = NA
      C1time = NA
      C2time = NA
      
      S = double(tpoints)
      L = double(tpoints)
      A1 = double(tpoints)
      A2 = double(tpoints)
      Y = double(tpoints)
      Ym = double(tpoints)
      pA1 = double(tpoints)
      pA2 = double(tpoints)
      pS = double(tpoints)
      chk = double(tpoints)
      ty = double(tpoints)
      Cense = double(tpoints)
      Cens  = double(tpoints)
      pcens = double(tpoints)
      pcense = double(tpoints)
      
      time = 1.0
      j = 1
      Ym[1]=0
      
      logitS = beta0 + beta1 * IT0
      pS[1] = 1/(1+(1/exp(logitS)))
      if (pS[1]==1) {S[1] = 1} else {S[1] = rbinom(1,1,pS[1])}
      L[1] = beta0 + beta1 * (1/log(T0))
      
      logitA1 = alphaA1_0 + alphaA1_2 * L[1]+ alphaA1_22 * S[1]
      pA1[1] = 1/(1+(1/exp(logitA1)))
      A1[1] = rbinom(1,1,pA1[1])
      
      logitA2 = alphaA2_0 + alphaA2_2 * L[1]+ alphaA2_22 * S[1]
      pA2[1] = 1/(1+(1/exp(logitA2)))
      A2[1] = rbinom(1,1,pA2[1])
      
      # generate Y's
      chk[1] = exp(psiA1*A1[1]+psiA2*A2[1])
      if (T0 > chk[1]) {
        Y[1]=0
        ty[1]=NA
      } else {
        Y[1]=1
        ty[1] = T0*exp(-psiA1*A1[1]-psiA2*A2[1])
      }
      
      j = 1
      Cense[1] = 0
      Cens[1] = 0
      pcens[1] = 0
      recurrent_A1_count <- 0
      recurrent_A2_count <- 0
      for (j in 2:tpoints){
        time = j
        if (!is.na(C2time) & j == C2time) {j = tpoints}
        if (!is.na(C1time) & j == C1time) {j = tpoints}
        Ym[j]=Y[j-1]
        
        logitS = beta0 + beta1 * IT0 + beta2 * A1[j-1] + beta3 * S[j-1]
        pS[j]=exp(logitS)/(1+exp(logitS))
        L[j] = beta0 + beta1 * (1/log(T0))+ beta2 * A1[j-1] + beta3 * L[j-1]
        
        if (pS[j] == 1) {S[j] = 1} else {S[j] = rbinom(1,1,pS[j])}
        if (A1[1] == 1) {
          A1_duration <- rtpois(1,5)
          A1[1:(1+A1_duration)] <- 1
          recurrent_A1_count <- recurrent_A1_count + 1
        }
        
        if (A1[j-1] == 0  && A1[1] != 1){
          logitA1 = alphaA1_0 + alphaA1_1 * A1[j-1] + alphaA1_2 * L[j] + alphaA1_3 * L[j-1] + alphaA1_22 * S[j] + alphaA1_32 * S[j-1] - recurrent_A1_count * 0.5
          if (logitA1 > 15){
            logitA1 <- 10
          }
          pA1[j]=exp(logitA1)/(1+exp(logitA1))
          pp= pA1[j]
          A1[j] = rbinom(1,1,pA1[j])
          if (A1[j] == 1) {
            A1_duration <- rtpois(1,5)
            A1[j:(j+A1_duration)] <- 1
            recurrent_A1_count <- recurrent_A1_count + 1
          }
        }
        if (recurrent_A1_count == 4){
          A1[j:tpoints] <- 0
        }
        
        if (A1[j]==1) {time=j}
        if (Ym[j]==1) {
          A1[j]=0
          S[j]=0
        }
        
        if (A2[1] == 1) {
          A2_duration <- rtpois(1,6)
          A2[1:(1+A2_duration)] <- 1
          recurrent_A2_count <- recurrent_A2_count + 1
        }
        if (A2[j-1] == 0 && A2[1] != 1){ 
          logitA2 = alphaA2_0 + alphaA2_1 * A2[j-1] + alphaA2_2 * L[j] + alphaA2_3 * L[j-1] + alphaA2_22 * S[j] + alphaA2_32 * S[j-1] - recurrent_A2_count * 0.5
          if (logitA2 > 15){
            logitA2 <- 10
          }
          pA2[j]=exp(logitA2)/(1+exp(logitA2))
          pp= pA2[j]
          A2[j] = rbinom(1,1,pA2[j])
          if (A2[j] == 1) {
            A2_duration <- rtpois(1,6)
            A2[j:(j+A2_duration)] <- 1
            recurrent_A2_count <- recurrent_A2_count + 1
          }
        }
        # if (A2[j-1] == 1){
        #   A2[j] <- 1
        # }
        if (recurrent_A2_count == 4){
          A2[j:tpoints] <- 0
        }
        if (A2[j]==1) {time=j}
        if (Ym[j]==1) {
          A2[j]=0
          S[j]=0
        }
        
        # generate Y's
        chk[j] = chk[j-1] + exp(psiA1*A1[j] + psiA2*A2[j])
        
        if (T0 > chk[j]) {
          Y[j]=0
          ty[j]=NA
        } else {
          Y[j]=1
          if (Ym[j]==1) {ty[j]<-ty[j-1]} else {ty[j]<-(j-1)+((T0-chk[j-1])*exp(-psiA1*A1[j]-psiA2*A2[j]))}
        }
        
        Cense[j] = Cense[1]
        Cens[j] = 0
        pcens[j] = 0
        pcense[j] = pcense[1]
        cens_1 = Cens[j]
        
      } #end j = 2 to tpoints
      
      if (is.na(ty[tpoints])) {
        C1time = tpoints
        C2time = tpoints
      } else {
        C1time = ceiling(ty[tpoints])
        C2time = ceiling(ty[tpoints])
      }
      
      time_used = min(C1time,C2time);
      
      
      output <- cbind(id, time = 0:(tpoints-1)/2, Y = Ym, A1 = A1[1:tpoints], A2 = A2[1:tpoints], L1 = L, L2 = S)
      
      outputx <- rbind(outputx, output)
      
      setTxtProgressBar(pb, id)
    } # end of id loop
    # outputx
    close(pb)
    dataset2 <- as.data.frame(outputx)
    write.csv(dataset2, file = paste("covid_simulation_data_", select.seed.val, ".csv", sep = ""), row.names = F)
    return(dataset2)
    cat("Output created for iteration", select.seed.val, "\n")
  }
  cat("Datasets being saved in:", getwd(), "\n")
  return(NULL)
}

