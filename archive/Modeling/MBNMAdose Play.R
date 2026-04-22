library(MBNMAdose)
library(ggdist)

# Load the SSRI dataset
data(ssri)

# Create a network object from the SSRI data
network <- mbnma.network(ssri)

# Print network summary to understand the data structure
print(network)
summary(network)

# Run the MBNMA model with the SSRI data

# Fits a Bayesian dose-response for model-based network meta-analysis (MBNMA) 
# that can account for multiple doses of different agents by applying a desired 
# dose-response function. Follows the methods of Mawdsley et al. (2016)

result <- mbnma.run(
  network = network,
  fun = dpoly(degree = 1),          
# Linear dose-response, can also use dlin(), dpoly(degree), dbs(), demax(), duser(), dmulti(), etc.
 
  method = "random",                
# Can also use "random" 
  
  regress = NULL,                   
# A formula of effect modifiers e.g. ~ Population + Age
  
  regress.effect = "common",        
# Can also "random" (assumed to be exchangeable versus Placebo throughout the network), 
# "agent" (assumed to be equal versus Placebo within each agent), 
# or "class" (assumed to be equal versus Placebo within each class)
  
  class.effect = list(),            
# e.g., list("dexamethasone" = "random"/"common")
  
  UME = FALSE,                      
# A boolean object to indicate whether to fit an Unrelated Mean Effects model 
# that does not assume consistency and so can be used to test if the consistency assumption is valid
  
  sdscale = FALSE,                  
# Logical object to indicate whether to write a model that specifies a reference SD 
# for standardising when modelling using Standardised Mean Differences
  
  cor = FALSE,                      
# Whether correlation should be modelled between relative effect dose-response parameters
  
  omega = NULL,                     
# A scale matrix for the inverse-Wishart prior for the covariance matrix used to model the 
                                    
# correlation between dose-response parameters. If left as NULL (the default) a diagonal matrix with
# elements equal to 100 is used.
  
  parameters.to.save = NULL,        
# A character vector containing names of parameters to monitor in JAGS
  
  pD = TRUE,                        
# If TRUE (the default) then adds the computation of pD, using the method
# of (Plummer 2008). If FALSE then uses the approximation of pD=var(deviance)/ 2 (often referred to as pV).
  
  likelihood = NULL,                
# Can take either "binomial", "normal" or "poisson". If left as NULL the 
# likelihood will be inferred from the data.
  link = NULL,                      
# Can take any link function defined within JAGS (e.g. "logit", "log", "probit", "cloglog"), be
# assigned the value "identity" for an identity link function, or be assigned the
# value "smd" for modelling Standardised Mean Differences using an identity link
# function. If left as NULL the link function will be automatically assigned based
# on the likelihood.
  
  priors = NULL,                    
# Can include a named list of parameter values (without indices) and replacement prior distribution. 
# values given as strings using distributions as specified in JAGS syntax. 
  
  n.iter = 20000,                   
# Number of total iterations per chain (including burn in; default: 20000)
  
  n.chains = 3,                     
# Number of Markov chains (default: 3)

  n.burnin = floor(20000/2),        
# length of burn in, i.e. number of iterations to discard at the beginning. Default is ‘n.iter/2“, 
# that is, discarding the first half of the simulations. If n.burnin is 0, jags() will run 100 iterations for adaption.    
  
  n.thin = max(1, floor((20000 - 10000)/1000)), 
# Thinning rate (default: max(1, floor((n.iter - n.burnin)/1000)))
  
  autojags = FALSE,                 
# A boolean value that indicates whether the model should be continually updated
#  until it has converged below a specific cutoff of Rhat
  
  Rhat = 1.05,                      
# A cutoff value for the Gelman-Rubin convergence diagnostic. Default is 1.05.

  n.update = 10,                    
# The maximum number of updates. Default is 10.

  model.file = NULL,

# The file path to a JAGS model (.jags file) that can be used to overwrite the JAGS
# model that is automatically written based on the specified options in MBNMAdose.
  
  jagsdata = NULL
# A named list of the data objects to be used in the JAGS model, if defining own
)

# Check model results
print(result)
summary(result)

# Plot dose-response curves
plot(result)

# Basic frequency table
table(df$variable)

# With proportions
prop.table(table(df$variable))

# Two-way table
table(ssri$agent)
