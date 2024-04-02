


//Data block: Define and name the size of each observed variable

data{    

   int N;             //Number of observations, i.e. trials
   int N_id;          //Number of individuals
   int id[N];         //Unique individual identification 
   int choice[N];        //Chosen option
  real reward[N];        //Amount of reward
   int blocktype[N];        //block
   int trialnum[N];        //trialnum

} 

//Parameter block: Define and name the size of each unobserved variable. 

parameters{ 
  
  //Learning parameters phi (updating rate) and lambda (inverse temperature)
  //parameters are on logit and log scale to incorporate individual random effects (see below)
  
   real logit_phi[1];   // One learning rate parameter      
   real log_L[1];           // One temperature parameter
   
  matrix[2, N_id] z_ID;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[2] sigma_ID;  //Standard deviation of learning parameters among individuals 
  cholesky_factor_corr[2] Rho_ID; //Cholesky factor for covariance of learning parameters among individuals
  
} 
   
transformed parameters{
  matrix[N_id, 2] v_ID;    
  v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)';    
}


model{

  matrix[N_id, 2] A; //Attraction matrix for 2 possible options (lakes)
  vector[N_id] lambda; //individual inv temp
  vector[N_id] phi; //individual learning rates
  //Assign priors
 
  logit_phi~  normal(0, 2); 
  
  
  log_L    ~  normal(0, 2);   
  
  //Define prior distribution of varying individual effects
  
  to_vector(z_ID) ~ normal(0, 1);   //This needs to be (0,1) prior because we want z-scores
  sigma_ID ~ exponential(1);         
  Rho_ID ~ lkj_corr_cholesky(4); 
  
  
  phi= rep_vector(0, N_id); 
  
  lambda= rep_vector(0, N_id); 
  
  //Initialize attraction scores
  for (i in 1:N_id) A[i, 1:2] = rep_vector(0.5, 2)'; 
  
  //Loop over Choices
  
  for (i in 1:N) {
  
  //Define and name local variables that update across choices
  
  //vector[4] pay;     //Payoffs    
  vector[2] p;       //Choice probabilities for 2 choices
 // real lambda;       //Inverse temp. on outcome scale   
  //real phi;          //Updating rate on outcome scale
 
  if (trialnum[i]==1){//reset prior in the begining of each block
  A[id[i],1:2]=rep_vector(0.5,2)';
  }
  
  //First, what is the log-probability of observed choice
  
  lambda[id[i]] = exp(log_L[1]  + v_ID[id[i], 1]);
  
  p = softmax(lambda[id[i]] * A[id[i], 1:2]' );      
  
  
  choice[i] ~ categorical(p);                                                                              
  
  //Second, update attractions conditional on observed choice
  
  phi[id[i]] = inv_logit(logit_phi[1] + v_ID[id[i], 2]);

  // Update attraction for all lakes

    A[ id[i] , choice[i] ] =  (1-phi[id[i]])*A[ id[i] , choice[i] ] + phi[id[i]]*reward[i];
    A[ id[i] , 3-choice[i] ] =  A[ id[i] ,3- choice[i] ] ;


  }
}

generated quantities{ 


  matrix[N_id, 2] A; //Attraction matrix for 2 possible options (lakes)
  vector[N_id] lambda; //individual inv temp
  vector[N_id] phi; //individual learning rates
  vector[N] log_lik; //log likelihood
  
  //Assign priors
 
  //logit_phi ~  normal(0, 2); 
  
  
 // log_L    ~  normal(0, 2);   
  
  //Define prior distribution of varying individual effects
  
//  to_vector(z_ID) ~ normal(0, 1);   //This needs to be (0,1) prior because we want z-scores
 // sigma_ID ~ exponential(1);         
//  Rho_ID ~ lkj_corr_cholesky(4); 
  
  phi= rep_vector(0, N_id); 
  lambda= rep_vector(0, N_id); 
  
  //Initialize attraction scores
  for (i in 1:N_id) A[i, 1:2] = rep_vector(0.5, 2)'; 
  
  //Loop over Choices
  
  for (i in 1:N) {
  
  //Define and name local variables that update across choices
  
 //vector[4] pay;     //Payoffs    
  vector[2] p;       //Choice probabilities for 2 choices
 //real lambda;       //Inverse temp. on outcome scale   
  //real phi;          //Updating rate on outcome scale
 
  if (trialnum[i]==1){//reset prior in the begining of each block
  A[id[i],1:2]=rep_vector(0.5,2)';
  }
  
  //First, what is the log-probability of observed choice
  
  lambda[id[i]] = exp(log_L[1]+v_ID[id[i], 1]);
  p = softmax(lambda[id[i]] * A[id[i], 1:2]' );      
  
  
  log_lik[i] = categorical_lpmf(choice[i] | p);                                                                              
  
  //Second, update attractions conditional on observed choice
  
  phi[id[i]] = inv_logit(logit_phi[1] + v_ID[id[i], 2]);

  // Update attraction for all lakes

    A[ id[i] , choice[i] ] =  (1-phi[id[i]])*A[ id[i] , choice[i] ] + phi[id[i]]*reward[i];
    A[ id[i] , 3-choice[i] ] =  A[ id[i] ,3- choice[i] ] ;


  }
}

