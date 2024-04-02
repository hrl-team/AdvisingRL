

RWmodelPE = function(data, par) {
  ll <- NA
  beta <- par[1] 
  alpha <- par[2]
 
  condition <-  as.numeric(as.character(data$blocktype))+1 # difficulty ???
  choice <- data$choice+1 #1 or 2 
  reward <- data$reward/10 # make 0.5 or 0 (in euro) to be 1 or 0 (yes or no)
  num_trials <-length(data$id) # 60trials
  x <-rep(0, num_trials)
  p_choice <- vector(mode = "numeric", num_trials) #initiate probability of choice with 0
 
  
  result <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("Q_choice", "Q_not_choice","p_choice","PE")
  
  colnames(result) <- x
  
  
  #initiate Q value matrix 2x4, with starting value of 0.5, to address each permutation of choice and condition (8)
  Q <- matrix(0.5, nrow = 2, ncol = 3, byrow = TRUE,
              dimnames = list(c("Choice 1", "Choice 2"),
                              c("difficulty 0 ", "difficulty 1", "difficulty 2")))
  PE<-0
  #for each rep, calculate the P choice (use exp function), PE and update the Q value
  for (RepNum in 1:num_trials) {
    
    # print(condition[RepNum])
    #use (3 - choice) to access the other choice in the same condition in the matrix:
    #3 - 1 = 2 or 3 - 2 = 1, gives the other option in the same condition
    p_choice[RepNum] <- exp(beta*Q[choice[RepNum],condition[RepNum]])/
      (exp(beta*Q[1,condition[RepNum]])+exp(beta*Q[2,condition[RepNum]]))
    
    result = rbind(result, data.frame(Q[choice[RepNum],condition[RepNum]],Q[3-choice[RepNum],condition[RepNum]],p_choice[RepNum],PE))
    
    #calculate prediction error
    PE <- reward[RepNum] - Q[choice[RepNum], condition[RepNum]]
    # x[RepNum] <- PE
  
    # Pe_df[RepNum][PE]<- PE
    #update Q value
   Q[choice[RepNum], condition[RepNum]] <- Q[choice[RepNum], condition[RepNum]] + alpha*PE
   
  
  }
  #minus maximum log likelihood, use sum and log functions
  # ll <- -sum(log(p_choice))
  
  #return ll
  return(result)
}

