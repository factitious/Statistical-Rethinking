# 3.1. Sampling from a grid-approximate posterior
# Compute the posterior for the globe tossing model, using grid approximation
# RCode 3.2

p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

# Sample from the posterior
samples <- sample( p_grid , prob=posterior, size=1e4 , replace=TRUE )
plot(samples)

# Compute a density estimate from theses samples
library(rethinking)
dens(samples)

#________________________________________________________________________________________________
#________________________________________________________________________________________________
# 3.2. Sampling to summarize
# Intervals of defined boundries

# the posterior probability that the proportion of water is less than 0.5. 
# Using the grid-approximate posterior, you can just add up all of the 
# probabilities, where the corresponding parameter value is less than 0.5

# RCode 3.6

# add up posterior probability where p < 0.5
sum( posterior[ p_grid < 0.8 ] )

# same calculation, using samples from the posterior:
# i.e. find the frequency of parameter values below 0.5
sum( samples < 0.5 ) / 1e4

# Using the same approach, you can ask how much posterior probability lies between 0.5 and 0.75:
sum( samples > 0.5 & samples < 0.75 ) / 1e4

#________________________________________________________________________________________________
# 3.2.2. Intervals of defined mass.
# i.e. Credible Interval (or Confidence Interval in FreqStats)
# These posterior intervals report two parameter values that contain between them a specified
# amount of posterior probability, a probability mass.

# Suppose for example you want to know the boundaries of the lower 80% posterior probability:
quantile( samples , 0.8 )

# Similarly, the middle 80% interval lies between the 10th percentile and the 90th percentile:
quantile( samples , c( 0.1 , 0.9 ) )

# Intervals of this sort, which assign equal probability mass to each tail, are very common
# in the scientific literature. We’ll call them percentile intervals (PI)
# These intervals do a good job of communicating the shape of a distribution, as long as the 
# distribution isn’t too asymmetrical. But in terms of supporting inferences about which parameters 
# are consistent with the data, they are not perfect.

# SEE BELOW + Book (pp. 56)

## 3W in 3Tosses + flat prior
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep(1,1000)
likelihood <- dbinom( 3 , size=3 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
samples <- sample( p_grid , size=1e4 , replace=TRUE , prob=posterior )

# Compute 50% percentile confidence interval
PI( samples , prob=0.5 )



HPDI( samples , prob=0.5 )

#________________________________________________________________________________________________
# 3.2.3. Point estimates. 

# The third and final common summary task for the posterior is to
# produce point estimates of some kind. Given the entire posterior distribution, what value
# should you report? This seems like an innocent question, but it is difficult to answer. The
# Bayesian parameter estimate is precisely the entire posterior distribution, which is not a single
# number, but instead a function that maps each unique parameter value onto a plausibility
# value. So really the most important thing to note is that you don’t have to choose a point estimate.
# It’s hardly ever necessary.

# See rest in BOOK
# Interesting yse if sapply for obtaining a list of loss values based on p_grid:
loss <- sapply( p_grid , function(d) sum( posterior*abs( d - p_grid ) ) )

#________________________________________________________________________________________________
#________________________________________________________________________________________________

# 3.3. Sampling to simulate prediction

# In this final section of the chapter, we’ll look at how to produce 
# simulated observations and how to perform some simple model checks.

#________________________________________________________________________________________________
# 3.3.1. Dummy data.
# RCode 3.20

dbinom( 0:2 , size=2 , prob=0.7 )

# Now we’re going to simulate observations, using these likelihoods. This is done by sampling
# from the distribution just described above. You could use sample to do this, but R
# provides convenient sampling functions for all the ordinary probability distributions, like
# the binomial. So a single dummy data observation of w can be sampled with:

rbinom( 1 , size=2 , prob=0.7 )
rbinom( 10 , size=2 , prob=0.7 )

dummy_w <- rbinom( 1e5 , size=2 , prob=0.7 )
table(dummy_w)/1e5

dummy_w <- rbinom( 1e5 , size=20 , prob=0.7 )
simplehist( dummy_w , xlab="dummy water count" )

#________________________________________________________________________________________________


w <- rbinom( 1e4 , size=9 , prob=0.6 )

# All you need to propagate parameter uncertainty into these predictions is replace the
# value 0.6 with samples from the posterior:

w <- rbinom( 1e4 , size=9 , prob=samples )

#________________________________________________________________________________________________
#________________________________________________________________________________________________

#3.5 Practice
#Easy
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
set.seed(100)
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )

#3E1. How much posterior probability lies below p = 0.2?
sum( samples < 0.2 ) / 1e4
# A = 5e-04 || 0.0005

#3E2. How much posterior probability lies above p = 0:8?
sum( samples > 0.8 ) / 1e4
# A = 0.1117

#3E3. How much posterior probability lies between p = 0:2 and p = 0:8?
sum( samples > 0.2 & samples < 0.8 ) / 1e4
# A = 0.8878

#3E4. 20% of the posterior probability lies below which value of p?
quantile(samples, 0.2)
# A = 0.52

#3E5. 20% of the posterior probability lies above which value of p?
quantile(samples, c(0.8, 1))
# A = 0.7567568 (-> 0.9769770)

#3E6. Which values of p contain the narrowest interval equal to 66% of the posterior probability?
HPDI( samples , prob=0.66 )
# A = 0.5205205 0.7847848 

#3E7. Which values of p contain 66% of the posterior probability, assuming equal posterior probability
# both below and above the interval?
PI(samples, prob = 0.66)
#A = 0.5005005 0.7687688 

#________________________________________________________________________________________________
# Medium

#3M1. Suppose the globe tossing data had turned out to be 8 water in 15 tosses. Construct the posterior
# distribution, using grid approximation. Use the same flat prior as before.
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 8 , size=15 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
set.seed(100)

#3M2. Draw 10,000 samples from the grid approximation from above. 
# Use the samples to calculate the 90% HPDI for p.
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )
HPDI(samples, 0.9)
#A = 0.3243243 -> 0.7157157 

#3M3. Construct a posterior predictive check for this model and data. This means simulate the distribution
# of samples, averaging over the posterior uncertainty in p. 
w <- rbinom( 1e4 , size=15 , prob=samples )
# What is the probability of observing 8 water in 15 tosses?
sort(table(w),decreasing=TRUE)/1e4 # / mean(w==8)
#A = 0.1484

# 3M4. Using the posterior distribution constructed from the new (8/15) data, now calculate the probability
# of observing 6 water in 9 tosses.
w_n <- rbinom(1e4, size = 9, prob = samples)
sort(table(w_n),decreasing=TRUE)/1e4
#A = 0.1824

#3M5. Start over at 3M1, but now use a prior that is zero below p = 0:5 and a constant above p = 0:5.
# This corresponds to prior information that a majority of the Earth’s surface is water. Repeat each
# problem above and compare the inferences. What difference does the better prior make? If it helps,
# compare inferences (using both priors) to the true value p = 0:7.

possible_p <- seq( from=0 , to=1 , length.out=1000 )
InfPrior <- rep( 1 , 1000 )
InfPrior <- ifelse(possible_p < .5, 0, 10)
likelihood_new <- dbinom( 8 , size=15 , prob=possible_p )
posterior_new <- likelihood_new * InfPrior
posterior_new <- posterior_new / sum(posterior_new)

#Draw 10,000 samples from the grid approximation from above. 
# Use the samples to calculate the 90% HPDI for p.
samples_new <- sample( possible_p , prob=posterior_new , size=1e4 , replace=TRUE )
HPDI(samples_new, 0.9)
# A =   0.5005005 -> 0.7067067 
# A_P = 0.3243243 -> 0.7157157 

# Construct a posterior predictive check for this model and data. This means simulate the distribution
# of samples, averaging over the posterior uncertainty in p. 
post_pred_distr <- rbinom( 1e4 , size=15 , prob=samples_new )
# What is the probability of observing 8 water in 15 tosses?
sort(table(post_pred_distr),decreasing=TRUE)/1e4 # / mean(w==8)
# A   = 0.162
# A_P = 0.1484

# Using the posterior distribution constructed from the new (8/15) data, now calculate the probability
# of observing 6 water in 9 tosses.
post_pred_distr2 <- rbinom(1e4, size = 9, prob = samples_new)
sort(table(post_pred_distr2),decreasing=TRUE)/1e4
# A = 0.2323
# A_P = 0.1824

#________________________________________________________________________________________________
# Hard

birth1 <- c(1,0,0,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,0,1,0,0,0,1,0,
            0,0,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,
            1,1,0,1,0,0,1,0,0,0,1,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,0,
            1,0,1,1,1,0,1,1,1,1)
birth2 <- c(0,1,0,1,0,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,0,
            1,1,1,0,1,1,1,0,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,1,1,
            0,0,0,1,1,1,0,0,0,0)

births <- data.frame(birth1,birth2)

#3H1. Using grid approximation, compute the posterior distribution for the probability of a birth
# being a boy. Assume a uniform prior probability. Which parameter value maximizes the posterior
# probability?

p_grid <- seq( from=0 , to=1 , length.out=1000 )
priorUniform <- rep(1,1000) 
lhBoy <- dbinom(round((sum(births$birth1)+sum(births$birth2))/2), size = length(births$birth1), prob = p_grid)
unstd.posteriorBoy <- lhBoy * priorUniform
posteriorBoy <- unstd.posteriorBoy/sum(unstd.posteriorBoy)
A_3H1 <- p_grid[ which.max(posteriorBoy) ]

#v2 - CORRECT

birthsTotal <- length(births$birth1)+length(births$birth2)
boysTotal <- sum(births)

lhBoy2 <- dbinom(boysTotal, birthsTotal, prob = p_grid)
unstd.posteriorBoy2 <- lhBoy2 * priorUniform
posteriorBoy2 <- unstd.posteriorBoy2/sum(unstd.posteriorBoy2)
A_3H1_v2 <- p_grid[ which.max(posteriorBoy2) ]


plot(posteriorBoy ~ p_grid, type = "l")
plot(posteriorBoy2 ~ p_grid, type = "l")

#3H2. Using the sample function, draw 10,000 random parameter values from the posterior distribution
# you calculated above. Use these samples to estimate the 50%, 89%, and 97% highest posterior
# density intervals.

samples <- sample( p_grid , prob=posteriorBoy2 , size=1e4 , replace=TRUE )

# 50% HPDI
HPDI(samples, 0.5)
# 0.5315315 0.5775776 

# 89% HPDI
HPDI(samples, 0.89)
# 5005005 0.6106106 

# 97% HPDI
HPDI(samples, 0.97)
# 0.4764765 0.6246246

#3H3. Use rbinom to simulate 10,000 replicates of 200 births. You should end up with 10,000 numbers,
# each one a count of boys out of 200 births. Compare the distribution of predicted numbers
# of boys to the actual count in the data (111 boys out of 200 births). There are many good ways to
# visualize the simulations, but the dens command (part of the rethinking package) is probably the
# easiest way in this case. Does it look like the model fits the data well? That is, does the distribution
# of predictions include the actual observation as a central, likely outcome?

post_pred_distr <- rbinom( 1e4 , size=200 , prob=samples )
sort(table(post_pred_distr),decreasing=TRUE)/length(post_pred_distr)
dens(post_pred_distr, adj = .1)
abline(v = boysTotal, col = "red")

#3H4. Now compare 10,000 counts of boys from 100 simulated first borns only to the number of boys
# in the first births, birth1. How does the model look in this light?

firstBoy <- sum(births$birth1)
ppd_firstBoy <- rbinom( 1e4 , size=100 , prob=samples )
sort(table(ppd_firstBoy),decreasing=TRUE)/length(ppd_firstBoy)
dens(ppd_firstBoy, adj = .1)
abline(v = firstBoy, col = "red")

#3H5. The model assumes that sex of first and second births are independent. To check this assumption,
# focus now on second births that followed female first borns. Compare 10,000 simulated counts
# of boys to only those second births that followed girls. To do this correctly, you need to count the
# number of first borns who were girls and simulate that many births, 10,000 times. Compare the
# counts of boys in your simulations to the actual observed count of boys following girls. How does the
# model look in this light? Any guesses what is going on in these data?

secondBirth_firstFem <- births$birth2[births$birth1==0] 
secondBoy_firstFem <- sum(secondBirth_firstFem)

# Direct way of doing this: 
secondBoy_firstFem <- sum(births$birth2[births$birth1==0]) 

ppd_secondBoy <- rbinom(1e4, size = length(secondBirth_firstFem), prob = samples)
sort(table(ppd_secondBoy),decreasing=TRUE)/length(ppd_secondBoy)
dens(ppd_secondBoy, adj = .1)
abline(v = secondBoy_firstFem, col = "red")

