####    Import Libraries    ####
library(insuranceData)
library(rcompanion)
library(tidyverse)
library(mixR)
library(actuar)
library(MASS)
library(ggplot2)
library(grid)
library(gridExtra)
library(egg)
library(stats)
library(boot)
library(dplyr)


#### Data Ingestion ####
data(dataCar)
data(AutoClaims)
p <- dataCar
c <- AutoClaims
claims_df <- c
policy_df <- p


#randomly sample from the claims in AutoClaims data and replace AGL in dataCar w. AGL in AutoDat

#set the seed here: this will ensure we get the same dataset every time we run this script 
#   Note: the dataset used in this analysis was created by sampling from the PAID column from table AutoClaims 
#     The seed is "reset" later on to ensure we don't get the same results every time we run the simulations. 
set.seed(123)
for( i in 1:nrow(policy_df) ){
  
  curr_num_claims <- policy_df[i,"numclaims"]
  if(curr_num_claims > 0){
    rand_row <- claims_df[ sample(nrow(claims_df), curr_num_claims, replace=TRUE), ]
    rand_row_sample <- rand_row$PAID
    policy_df[i, "claimcst0"] <- sum(rand_row_sample)
  }
}


policy_df_static <- policy_df
non_zero_claims <- subset(policy_df, policy_df$numclaims>0)
zero_claims <- subset(policy_df, policy_df$numclaims<=0)
final_df <- rbind(non_zero_claims, zero_claims)



####    Transform Data & Calculate MLE's for Lognormal Severity Distribution   ####
#split data into 2004 and 2005 Policies: first subset only rows with exposure = 365/365.25
fully_exp_df <- subset(final_df, final_df$exposure >= 0.999)
#Assune that 70% of the policies in the fully_exp_df created above are from 2005 (fully exposed policies)
sample_05 <- fully_exp_df[ sample(nrow(fully_exp_df), 0.70*nrow(fully_exp_df), replace=FALSE), ]
#get indices corresponding to rows found from sampling done above
indices <- sample(1:nrow(sample_05), nrow(sample_05))
#Take the DF w. all policies fully exposed (mixture of 04/05 policies) and remove the rows that correspond to 05 policies
df_04 <- fully_exp_df[-indices,]
#Repeat above process this time remove the rows corresponding to 04 policies 
indices_04 <- sample(1:nrow(df_04), nrow(df_04))
df_05 <- final_df[-indices_04,]


#filter 05 data so we only have policies with exposure >= 0.90     
df_05 <- subset(df_05, df_05$exposure>=0.90)
df_05_nonzero <- subset(df_05, df_05$clm == 1)
#Limit the claim loss amount to < $20K
df_05_nonzero <- subset(df_05_nonzero, df_05_nonzero$claimcst < 1.7e4)


#calculate average claims count per policy for policies effective in the year 2005
avg_clm_count <- sum(df_05$numclaims)/nrow(df_05)
avg_clm_count_per_pol <- sum(df_05$numclaims)/nrow(df_05)
#non_zero_claims_05 <- subset(df_05, df_05$numclaims > 0)
total_pols_05 <- nrow(df_05)
#non_zero_claims_df <- subset(df_05, claimcst0 > 0)
mle_mu <- mean(log(df_05_nonzero$claimcst0))
mle_sigma <- sd(log(df_05_nonzero$claimcst0))


####  Poisson Process Example ####
#this is lambda
rate <- 0.50
#this is on the interval (0,t]
t <- 30
#draw from a poisson w. parameter lambda*t
poiss_draw <- rpois(1, rate*t)
#draw from random uniform on interval (0,t]
unifs <- runif(poiss_draw, 0, t)
#sort draws from uniform above
arrivals <- sort(unifs)
df <- as.data.frame(arrivals)
ggplot(data = df) + geom_step(aes(c(1:poiss_draw), df$arrivals)) + labs(title= "Poisson Process Example", x="Arrival Times", y="Number of Arrivals" )+ geom_point(aes(c(1:poiss_draw), df$arrivals))+ theme(plot.title = element_text(hjust = 0.5))


#this essentially "resets" the seed by taking the current time-stamp (changes every second)
#  this is needed b.c. we want to set the seed above when we create/transform our original data HOWEVER
#   We want the seed (initiliazation of the psuedo-random number generator) to be random when we run our simulations 
#   Note: there is another copy of this script that DOES NOT "reset" the seed - therefore you get the samwe results every time. 
set.seed(Sys.time())

####    EDA   ####
#plot hisogram of Exposures
ggplot(final_df) + geom_histogram(aes(x=exposure), bins=75) + labs(title= "Exposure: Auto Policies 2004/2005")+ theme(plot.title = element_text(hjust = 0.5)) 


ggplot(df_05_nonzero) + geom_histogram(aes(x=claimcst0), bins=400) + labs(title= "Losses: Auto Policies 2005")+ theme(plot.title = element_text(hjust = 0.5))  -> p1

ggplot(df_05_nonzero) + geom_density(aes(x=log(claimcst0)),color="black", fill="lightblue") + labs(title= "Log(Losses): Auto Policies 2005",xlab=c("Log Claim Loss Size"))+ theme(plot.title = element_text(hjust = 0.5))  -> p2

log_loss <- as.data.frame(log(df_05_nonzero$claimcst0))
colnames(log_loss) <- c("log_losses")

ggplot(data=log_loss, aes(sample = log_losses)) + stat_qq() + stat_qq_line() + labs(title= "Log(Losses) Q-Q Plot: Auto Policies 2005",xlab=c("Log Claim Loss Size"))+ theme(plot.title = element_text(hjust = 0.5)) -> p3

ggarrange(p1,p2,p3 ,nrow=3)




####    Claim Simulation - Severity Lognormal   ####
claim_simulation <- function(avg_claim_count, total_policies, mean_sev_param, sd_sev_param){
  pois_rate <- avg_claim_count
  t <- total_policies
  pois_draw <- rpois(1, pois_rate*t)
  sim_agl_vals <- rlnorm(pois_draw, meanlog = mean_sev_param, sdlog = sd_sev_param)
  sim_agl <- sum(sim_agl_vals)
  
  return(sim_agl)
}

simulation_output <- replicate(1e4,claim_simulation(avg_clm_count,total_pols_05,mle_mu,mle_sigma))

#plotNormalHistogram(simulation_output, main="2005: Simulated Aggregate Losses $'s")
hist(simulation_output, main="2005: Simulated Aggregate Losses $'s", breaks=100)



####  Compare Simulation I to Actuals   ####
simulated_annual_losses <- mean(simulation_output)
actual_annual_losses <- sum(df_05$claimcst0)
delta <- actual_annual_losses - simulated_annual_losses
print(paste("Expected Annual Losses = ",round(simulated_annual_losses,digits=0)))
print(paste("Actual Annual Losses = ", round(actual_annual_losses, digits=0)))
print(paste("Pct Error = ",round(100*delta/simulated_annual_losses, digits=2)))




####  Bootstrap Likelihood Ratio Test to Find Number of Components in Finite Mixture Model    ####

##### Bootstrap Likelihood Ratio Test for Finite Mixture Models (confirm # components)
#calculate average claims count per policy for policies effective in the year 2005
#log_agl <- log(df_05_nonzero$claimcst0)
#find number of components (either one or two in this example)
# agl_output1 <- bs.test(log(df_05_nonzero$claimcst0), ncomp = c(1, 2), family = c("normal", "weibull", "gamma",
     #                                                                             "lnorm"), B = 100, ev = FALSE, mstep.method = c("bisection", "newton"),
    #                    init.method = c("kmeans", "hclust"), tol = 1e-06, max_iter = 500)


# #REJECT null - number of components is > 1
 #print(paste("P-value for 1-2 Components = ", agl_output1[1] ))

 #find number of components (either two or three in this example)
# agl_output2 <- bs.test(log(df_05_nonzero$claimcst0), ncomp = c(2, 3), family = c("normal", "weibull", "gamma",
     #                                                                             "lnorm"), B = 100, ev = FALSE, mstep.method = c("bisection", "newton"),
  #                      init.method = c("kmeans", "hclust"), tol = 1e-06, max_iter = 500)

# #REJECT null hypothesis - this implies > 2 components
 #print(paste("P-value for 2-3 Components = ", agl_output2[1] ))

 #find number of components (either three or four in this example)
# agl_output3 <- bs.test(log(df_05_nonzero$claimcst0), ncomp = c(3, 4), family = c("normal", "weibull", "gamma",
        #                                                                          "lnorm"), B = 100, ev = FALSE, mstep.method = c("bisection", "newton"),
  #                      init.method = c("kmeans", "hclust"), tol = 1e-06, max_iter = 500)

# #FAIL TO REJECT the Null Hypothesis - this implies we can assume a mixture with 3 components
 #print(paste("P-value for 3-4 Components = ", agl_output3[1] ))
#





#### fit mixture model of Lognormals with 3 components found using the Bootstrap Likelihood Ratio Test above  ####
log_agl_05 <- log(df_05_nonzero$claimcst0)
fit_agl_mix3 <- mixfit(log_agl_05, ncomp=3, family="lnorm") 
plot(fit_agl_mix3, main="Auto-Claims Log(Losses) Mixture Model: Three Component LogNormal Mixture")



####   Simualtion II: Finite Mixture Model    ####
#respective weights for each of the three components to our mixture model
# weight_c1 <- fit_agl_mix3$pi[1]
# weight_c2 <- fit_agl_mix3$pi[2]
# weight_c3 <- fit_agl_mix3$pi[3]
# 
# #respective means for each of the three components to our mixture model
# mu_c1 <- fit_agl_mix3$mu[1]
# mu_c2 <- fit_agl_mix3$mu[2]
# mu_c3 <- fit_agl_mix3$mu[3]
# 
# #respective means for each of the three components to our mixture model
# sd_c1 <- fit_agl_mix3$sd[1]
# sd_c2 <- fit_agl_mix3$sd[2]
# sd_c3 <- fit_agl_mix3$sd[3]

# #Create Compound Poisson Distribution using package
# severity_mixture_closed <- rcomppois(1e5, avg_claim_count_per_pol*total_pols_05,
#                                      rmixture(probs = fit_agl_mix3$pi,
#                                               expression(rlnorm(fit_agl_mix3$mu[1], fit_agl_mix3$sd[1]),
#                                                          rlnorm(fit_agl_mix3$mu[2], fit_agl_mix3$sd[2]),
#                                                          rlnorm(fit_agl_mix3$mu[3], fit_agl_mix3$sd[3]))))




claim_simulation_mix <- function(avg_claim_count, total_policies){
  pois_rate <- avg_claim_count
  t <- total_policies
  pois_draw <- rpois(1, pois_rate*t)
  sim_mix_agl_vals <- rmixlnorm(pois_draw, fit_agl_mix3$pi, fit_agl_mix3$mu, fit_agl_mix3$sd)
  exp_sim_mix_agl_vals <- exp(sim_mix_agl_vals)
  sim_agl_sum_mix <- sum(exp_sim_mix_agl_vals)
  
  return(sim_agl_sum_mix)
}

simulation_output_mix <- replicate(1e4, claim_simulation_mix(avg_clm_count_per_pol, total_pols_05))
#plotNormalHistogram(simulation_output_mix, main="2005: Simulated Aggregate Losses (Mixture Model)")
min_sim_mix <- min(simulation_output_mix)
hist(simulation_output_mix, breaks = 100, main="2005: Simulated Aggregate Losses (Mixture Model)", xlim=range(min_sim_mix:2.5e6) )




####      Compare Results From Mixrture Model Simulation to Actual    ####
simulated_annual_losses_mix <- mean(simulation_output_mix)
actual_annual_losses <- sum(df_05$claimcst0)
delta <- actual_annual_losses - simulated_annual_losses_mix
print(paste("Expected Annual Losses = ",round(simulated_annual_losses_mix,digits=0)))
print(paste("Actual Annual Losses = ", round(actual_annual_losses, digits=0)))
print(paste("Pct Error = ",round(100*delta/simulated_annual_losses_mix, digits=2)))

