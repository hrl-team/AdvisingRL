library("tidyverse")
library("ggplot2")
library("RColorBrewer")
library("MetBrewer")
library("tidyverse")
library("rstan")
library("ggthemes")
RWmodelPE_SIM <- function(num_trials, alphas, betas, num_iterations) {
  num_difficulties <- 3
  trials_per_difficulty <- num_trials / num_difficulties
  
  mean_rewards <- matrix(c(6.5, 6, 5.5, 3.5, 4, 4.5), nrow = 2, ncol = 3, byrow = TRUE,
                         dimnames = list(c("Choice 1", "Choice 2"),
                                         c("difficulty 1", "difficulty 2", "difficulty 3")))
  sd_reward <- 1.7
  
  all_results <- data.frame()
  sim_id <- 0
  
  for (alpha in alphas) {
    for (beta in betas) {
      for (iteration in 1:num_iterations) {
        sim_id <- sim_id + 1
        condition <- rep(1:num_difficulties, each = trials_per_difficulty)
        trial_numbers <- rep(1:trials_per_difficulty, num_difficulties)
        
        choice <- rep(NA, num_trials)
        reward <- rep(NA, num_trials)
        
        p_choice <- vector(mode = "numeric", num_trials)
        
        result <- data.frame("trial" = numeric(num_trials),
                             "Q_choice" = numeric(num_trials),
                             "Q_not_choice" = numeric(num_trials),
                             "p_choice" = numeric(num_trials),
                             "PE" = numeric(num_trials),
                             "difficulty" = numeric(num_trials),
                             "actual_choice" = numeric(num_trials),
                             "alpha" = numeric(num_trials),
                             "beta" = numeric(num_trials),
                             "iteration" = numeric(num_trials),
                             "SimID" = numeric(num_trials),
                             "reward" = numeric(num_trials),
                             "choicePoS" = numeric(num_trials))
        
        Q <- matrix(5, nrow = 2, ncol = 3, byrow = TRUE, 
                    dimnames = list(c("Choice 1", "Choice 2"),
                                    c("difficulty 1", "difficulty 2", "difficulty 3")))
        PE <- 0
        
        for (RepNum in 1:num_trials) {
          # Stabilized Softmax Calculation
          max_Q <- max(Q[1, condition[RepNum]], Q[2, condition[RepNum]])
          exp_Q1 <- exp(beta * (Q[1, condition[RepNum]] - max_Q))
          exp_Q2 <- exp(beta * (Q[2, condition[RepNum]] - max_Q))
          prob_choice_1 <- exp_Q1 / (exp_Q1 + exp_Q2)
          
          # Simulate choice based on the computed probabilities
          choice[RepNum] <- sample(c(1, 2), 1, prob = c(prob_choice_1, 1 - prob_choice_1))
          p_choice[RepNum] <- ifelse(choice[RepNum] == 1, prob_choice_1, 1 - prob_choice_1)
          
          # Simulate the reward for the chosen option
          reward[RepNum] <- min(8, max(2, round(rnorm(1, mean_rewards[choice[RepNum], condition[RepNum]], sd_reward))))
          
          # Calculate the prediction error
          PE <- reward[RepNum] - Q[choice[RepNum], condition[RepNum]]
          
          # Update the Q-value for the chosen option
          Q[choice[RepNum], condition[RepNum]] <- Q[choice[RepNum], condition[RepNum]] + alpha * PE
          
          # Determine if the actual choice was the one with the higher mean reward
          actual_choice <- ifelse(mean_rewards[choice[RepNum], condition[RepNum]] > mean_rewards[3 - choice[RepNum], condition[RepNum]], 1, 0)
          
          # Store the results of the current trial
          result[RepNum, ] <- c(trial_numbers[RepNum],
                                Q[choice[RepNum], condition[RepNum]], 
                                Q[3 - choice[RepNum], condition[RepNum]], 
                                p_choice[RepNum], PE, condition[RepNum], actual_choice, alpha, beta, iteration, sim_id, reward[RepNum], choice[RepNum])
        }
        
        # Append the results of the current iteration to the comprehensive results data frame
        all_results <- rbind(all_results, result)
      }
    }
  }
  
  return(all_results)
}

# Generate:
num_trials <- 60
alphas <- rbeta(9, shape1 = 1.1, shape2 = 1.1)
betas <- rgamma(9, shape = 1.2, scale = 5)
num_iterations <- 1

result <- RWmodelPE_SIM(num_trials, alphas, betas, num_iterations)

# head(result, n = 60)

DtAvsI <- result %>% group_by(trial, difficulty) %>% summarise(Accuracy=mean(actual_choice),
                                                                        AccuracySD=sd(replicate(1000, mean(sample(actual_choice, replace=T))))
                                                                        ) %>% distinct()

DtAvsI$blocktype <- factor(DtAvsI$difficulty)

ggplot() +
           geom_errorbar(data=DtAvsI, aes(x=trial, y=Accuracy, ymax=Accuracy + AccuracySD, ymin=Accuracy - AccuracySD), color="black", width=0.3) +
           geom_line(data=DtAvsI, aes(x=trial, y=Accuracy), color="black", linewidth=3, size=2) +
           geom_point(data=DtAvsI, aes(x=trial, y=Accuracy), fill="black", color="black", size=3, stroke=3, shape=21) +

           geom_hline(yintercept=0.5, linetype=3, linewidth=2) +
          scale_fill_manual(aesthetics = c("color", "fill"), values = canva_palettes$`Primary colors with a vibrant twist`[1:3]) +
  facet_wrap(vars(blocktype)) +
  theme_bw() +
  theme(legend.text = element_text(angle=0, face="plain", colour="black", size=30, family="arial"),
        legend.position.inside = 0,
        legend.key.size = unit(4, "line"),
        strip.text.x =  element_text(angle=0, face="plain", colour="black", size="38", family="arial"),
        strip.background = element_blank(),
        panel.spacing = unit(4, "lines"), #space between plots
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line=element_line(size = 1), 
        panel.border=element_blank(), 
        axis.text.x=element_text(size=28, family="arial"), axis.title.x=element_text(size=38, family="arial"),
        axis.text.y=element_text(size=28, family="arial"), axis.title.y=element_text(size=38, family="arial"),
        plot.title=element_text(size=22, family="arial", face="plain", color="black")) +
labs(face="plain", family="arial",
     x = "Trial", 
     y = "P(correct)")


#recover parameters from the sim
####RL: invoke STAN
##### 
#rename to match original script its faster

result$blocktype<-result$difficulty
result$id<-result$SimID
result$choice<-result$actual_choice
result$trialnum<-result$trial


dat <- as.list(result[,c("id", "choice","reward", "blocktype","trialnum")])
K=1
for (sub_id in unique(dat$id)){
  dat$id_2[dat$id==sub_id]=K
  K=K+1
}
dat$id=dat$id_2
dat$reward=dat$reward/10
dat$choice=dat$choice+1
dat$N <- nrow(result)
dat$N_id <- length(unique(result$id))

sm <- stan_model(file = "/home/subco/NEW_WORK/UriAdvice/RW_RL.stan") #MAKE SURE YOU SET THE CORRECT PATH TO THE STAN FILE
m <- sampling(sm, data = dat, iter = 4000, cores = 4, chains = 4, refresh = 1, control = list(adapt_delta = 0.9, max_treedepth = 14))

LR_summary <- summary(m, pars = c( "phi"), probs = c(0.1, 0.9))$summary
Temp_summary <- summary(m, pars = c("lambda"), probs = c(0.1, 0.9))$summary

RF_Param<-data.frame(LR=LR_summary[,1],Temp=Temp_summary[,1],id=seq(1,K-1))

#Check if recovery is good

result <- merge(result, RF_Param, by = "id")

ggplot(result, aes(x=scale(alpha), y=scale(LR))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Alpha", y = "Recovered Alpha") +
  theme_minimal()

  ggplot(result, aes(x=scale(beta), y=scale(Temp))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Beta", y = "Recovered Beta") +
  theme_minimal()

  cor.test(result$alpha, result$LR)
  cor.test(result$beta, result$Temp)

#check the quality of recovery also with ICC
library(psych)
ICC(data.frame(scale(result$alpha), scale(result$LR)))
ICC(data.frame(scale(result$beta), scale(result$Temp)))
