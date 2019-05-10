# Gaussian distribution by addition
pos <- replicate( 1e7 , sum( runif(16,-1,1) ) )
plot(density(pos))

# Gaussian distribution by multiplication ()
prod( 1 + runif(12,0,0.1) )

growth <- replicate( 1000000 , prod( 1 + runif(12,0,0.1) ) )
dens( growth , norm.comp=TRUE )

big <- replicate( 10000 , prod( 1 + runif(12,0,0.5) ) )
small <- replicate( 10000 , prod( 1 + runif(12,0,0.01) ) )
dens(big, norm.comp = T)
dens(small, norm.comp = T)
# Small effects that multiply together tend are approximately additive,
# and so they tend to stabilize on Gaussian distributions.

# Gaussian distribution by log-multiplication
log.big <- replicate( 100000 , log(prod(1 + runif(12,0,0.5))) )
dens(log.big, norm.comp = T)

# From model definition to Bayes’ theorem.
w <- 6; n <- 9;
p_grid <- seq(from=0,to=1,length.out=1000)
likelihood <- dbinom(w,n,p_grid)
prior <- dunif(p_grid,0,1)
posterior <- likelihood * prior
posterior <- posterior/sum(posterior)
plot(posterior)


#4.3. A Gaussian model of height

# Load data
# Partial census data for the Dobe area !Kung San, compiled from 
# interviews conducted by Nancy Howell in the late 1960s

data(Howell1) 
d <- Howell1

# Filter data by age (i.e. adults only)
d2 <- d[d$age >= 18,]

# Plot the distribution of height
dens(d2$height)

# Plot the prior for the mean
curve( dnorm( x , 178 , 20 ) , from=100 , to=250 )

# Plot the prior for the standard deviation 
# (The std prior is a truly flat prior, a uniform one, that functions
# just to constrain  to have positive probability between zero and 50cm)
curve( dunif( x , 0 , 50 ) , from=-10 , to=60 )

# Prior probability distribution of height

# You didn’t specify a prior probability distribution of heights
# directly, but once you’ve chosen priors for mu (mean) and sigma (std), these 
# imply a prior distribution of individual heights. 
# You can quickly simulate heights by sampling from the prior, like you would 
# from the posterior back in Chapter 3. Remember, every posterior is also potentially a
# prior for a subsequent analysis, so you can process priors just like posteriors.

sample_mu <- rnorm( 1e4 , 178 , 20 )
sample_sigma <- runif( 1e4 , 0 , 50 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma )
dens( prior_h )

# 4.3.3. Grid approximation of the posterior distribution. (pp. 84)

mu.list <- seq( from=140, to=160 , length.out=200 )
sigma.list <- seq( from=4 , to=9 , length.out=200 )
post <- expand.grid( mu=mu.list , sigma=sigma.list )
post$LL <- sapply( 1:nrow(post) , function(i) sum( dnorm( #LL - log-likelihood. 
  # Everything needs to be done on the log scale. Otherwise rounding error will quickly make all of the posterior probabilities zero
  # Sapply passes the unique combination of [mu] and [sigma] on each row of post to a function that computes 
  # the log-likelihood of each observed height, and adds all of these log-likelihoods together (sum)
  d2$height ,
  mean=post$mu[i] ,
  sd=post$sigma[i] ,
  log=TRUE ) ) )
post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) +
  dunif( post$sigma , 0 , 50 , TRUE )
post$prob <- exp( post$prod - max(post$prod) )

# Contour plot (of the posterior distribution)
contour_xyz( post$mu , post$sigma , post$prob )

# Heat map (of the posterior distribution)
image_xyz( post$mu , post$sigma , post$prob )

#4.3.4. Sampling from the posterior.
# The only new trick is that since there are two parameters, and
# we want to sample combinations of them

# randomly sample row numbers in post in proportion to the values in post$prob
sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       prob=post$prob ) 

# pull out the parameter values on those randomly sampled rows
sample.mu <- post$mu[ sample.rows ]
sample.sigma <- post$sigma[ sample.rows ]

plot( sample.mu , sample.sigma , cex=1 , pch=15 , col=col.alpha(rangi2,0.3) )

# Marginal posterior densities
# marginal = “averaging over the other parameters.”
dens( sample.mu )
dens( sample.sigma)

# Highest Posterior Density Intervals
HPDI (sample.mu)
HPDI (sample.sigma)

# Overthinking: Sample size and the normality of sigma’s posterior.

#sample 20 random heights
d3 <- sample( d2$height , size=20 )

mu.list <- seq( from=150, to=170 , length.out=200 )
sigma.list <- seq( from=4 , to=20 , length.out=200 )
post2 <- expand.grid( mu=mu.list , sigma=sigma.list )
post2$LL <- sapply( 1:nrow(post2) , function(i)
  sum( dnorm( d3 , mean=post2$mu[i] , sd=post2$sigma[i] ,
              log=TRUE ) ) )
post2$prod <- post2$LL + dnorm( post2$mu , 178 , 20 , TRUE ) +
  dunif( post2$sigma , 0 , 50 , TRUE )
post2$prob <- exp( post2$prod - max(post2$prod) )
sample2.rows <- sample( 1:nrow(post2) , size=1e4 , replace=TRUE ,
                        prob=post2$prob )
sample2.mu <- post2$mu[ sample2.rows ]
sample2.sigma <- post2$sigma[ sample2.rows ]
plot( sample2.mu , sample2.sigma , cex=0.5 ,
      col=col.alpha(rangi2,0.1) ,
      xlab="mu" , ylab="sigma" , pch=16 )

# marginal posterior density for mu, averaging over sigma:
dens( sample2.sigma , norm.comp=TRUE )

#4.3.5. Fitting the model with map.(RCode 4.24 and below)

flist <- alist(
                height ~ dnorm( mu , sigma ) ,
                mu ~ dnorm( 178 , 20 ) ,
                sigma ~ dunif( 0 , 50 )
)

m4.1 <- map( flist , data=d2)
precis(m4.1)

# These numbers provide Gaussian approximations for each parameter’s marginal distribution.
# This means the plausibility of each value of mu, after averaging over the plausibilities of each
# value of sigma, is given by a Gaussian distribution with mean 154.6 and standard deviation 0.4.

# Building a more informative prior into map
m4.2 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu ~ dnorm( 178 , 0.1 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d2 )
precis( m4.2 )

# Just like a mean and standard deviation (or its square, a variance) are sufficient to
# describe a one-dimensional Gaussian distribution, a list of means and a matrix of variances
# and covariances are sufficient to describe a multi-dimensional Gaussian distribution.
vcov( m4.1 )

# A variance-covariance matrix can be factored into two elements"
#1. A vector of variances for the parameters
diag( vcov( m4.1 ) )

#2. A correlation matrix that tells us how changes in any parameter lead to correlated changes in the others 
cov2cor( vcov( m4.1 ) )

# To get samples from the posterior:
# Instead of sampling single values from a simple Gaussian distribution, we sample 
# vectors of values from a multi-dimensional Gaussian distribution.

post <- extract.samples( m4.1 , n=1e4 )
head(post)

plot(post)
dens(post)

#4.4 Adding a predictor
plot( d2$height ~ d2$weight )

# fit model (now with predictor variable)
m4.3 <- map(
  alist(
    height ~ dnorm( mu , sigma ) , # height ~ dnorm( a + b*weight , sigma ) - embeds the linear model into the LH
    mu <- a + b*weight ,
    a ~ dnorm( 178 , 100 ) ,
    b ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d2 )


# The parameter mu is no longer really a parameter here, 
# because it has been replaced by the linear model, a+b*weight.
# So there is a prior for the parameter a now, but not one for mu, 
# since mu is defined by the linear model instead

# Starting values for MAP (determines where the function starts 'climbing' the posterior)
# You will need to run this before you fit the model (as map requires 'start')
# Does not appear to affect precis of model. Why?

start <- list(
  a=mean(d2$height),
  b=0,
  sigma=sd(d2$height)
)

precis(m4.3)


# The numbers in the default precis output aren’t sufficient to describe  the quadratic 
# posterior completely. For that, we also require the variance-covariance matrix.
precis( m4.3 , corr=TRUE )

# In more complex models, strong correlations like this can make it difficult to 
# fit the model to the data. So we’ll want to use some golem engineering tricks to avoid it, when possible.

# CENTERING:
# Centering is the procedure of subtracting the mean of a variable from each value

# To create a centered version of the weight variable:
d2$weight.c <- d2$weight - mean(d2$weight)

# Replace d2$weight with d2$wight.c in the model fit:
m4.4 <- map(
             alist(
               height ~ dnorm( mu , sigma ) ,
               mu <- a + b*weight.c ,
               a ~ dnorm( 178 , 100 ) ,
               b ~ dnorm( 0 , 10 ) ,
               sigma ~ dunif( 0 , 50 )
             ) ,
             data=d2 )

precis( m4.4 , corr=TRUE )

# The estimates for b and sigma are unchanged (within rounding error), but the estimate for a is
# now the same as the average height value in the raw data. Try it yourself: mean(d2$height).
# And the correlations among parameters are now all zero. What has happened here?

# The estimate for the intercept, a, still means the same thing it did before: the expected
# value of the outcome variable, when the predictor variable is equal to zero. But now the
# mean value of the predictor is also zero. So the intercept also means: the expected value of
# the outcome, when the predictor is at its average value. This makes interpreting the intercept
# a lot easier

#4.4.3.2 Plotting posterior inference against the data.

# To superimpose the MAP values for mean height over the actual data:
plot( height ~ weight , data=d2 )
abline( a=coef(m4.3)["a"] , b=coef(m4.3)["b"] )

#4.4.3.3. Adding uncertainty around the mean.

# To better appreciate how the posterior distribution contains lines, extract some samples from the model:
post <- extract.samples( m4.3 )

# The following code extracts the first 10 cases and re-estimates the model:
N <- 1000
dN <- d2[ 1:N , ]
mN <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*weight ,
    a ~ dnorm( 178 , 100 ) ,
    b ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) , data=dN )

# Now let’s plot 20 of these lines, to see what the uncertainty looks like.

# extract 20 samples from the posterior
post <- extract.samples( mN , n=20 )
# display raw data and sample size
plot( dN$weight , dN$height ,
      xlim=range(d2$weight) , ylim=range(d2$height) ,
      col=rangi2 , xlab="weight" , ylab="height" )
mtext(concat("N = ",N))
# plot the lines, with transparency
for ( i in 1:20 )
  abline( a=post$a[i] , b=post$b[i] , col=col.alpha("black",0.3) )

#4.4.3.4. Plotting regression intervals and contours.

mu_at_50 <- post$a + post$b * 50

dens( mu_at_50 , col=rangi2 , lwd=2 , xlab="mu|weight=50" )
# The quadratic approximate posterior distribution of the mean height, mu, when
# weight is 50 kg. This distribution represents the relative plausibility of 
# different values of the mean.

HPDI( mu_at_50 , prob=0.89 )

# That’s good so far, but we need to repeat the above calculation for every weight value on
# the horizontal axis, not just when it is 50 kg. We want to draw 89% HPDIs around the MAP slope

mu <- link( m4.3 )
str(mu)

# What link will do is take your map model fit, sample from the posterior distribution,
# and then compute mu for each case in the data and sample from the posterior distribution.
# Each row is a sample from the posterior distribution
# Each column is a case (row) in the data

# So above we have a distribution of mu for each individual in the original data. 
# We actually want something slightly different: a distribution of mu for each unique weight value 
# on the horizontal axis. 
# It’s only slightly harder to compute that, by just passing link some new data:

# Define sequence of weights to compute predictions for.
# These values will be on the horizontal axis

#weight.seq <- seq( from=25 , to=70 , by=1 ) 
# unique(m4.3@data[["weight"]]) - to include each unique weight value
weight.seq <- seq( round(min(d2$weight)) , round(max(d2$weight)) , by=1 ) 

# Use link to compute mu for each sample from posterior
# and for each weight in weight.seq
mu <- link( m4.3 , data=data.frame(weight=weight.seq) )
str(mu)


# use type="n" to hide raw data
plot( height ~ weight , d2 , type="n" )
# loop over samples and plot each mu value
for ( i in 1:100 )
  points( weight.seq , mu[i,] , pch=16 , col=col.alpha(rangi2,0.1) )

# At each weight value in weight.seq, a pile of computed mu values are shown. 
# Each of these piles is a Gaussian distribution

# The final step is to summarize the distribution for each weight value. We’ll use apply,
# which applies a function of your choice to a matrix.

# summarize the distribution of mu
mu.mean <- apply( mu , 2 , mean ) # compute the mean of each column (dimension “2”) of the matrix mu
mu.HPDI <- apply( mu , 2 , HPDI , prob=0.89 )

# plot raw data
# fading out points to make line and interval more visible
plot( height ~ weight , data=d2 , col=col.alpha(rangi2,0.5) )
# plot the MAP line, aka the mean mu for each weight
lines( weight.seq , mu.mean )
# plot a shaded region for 89% HPDI
shade( mu.HPDI , weight.seq )


# 4.4.3.5. Prediction intervals.
# What you’ve done so far is just use samples from the posterior to visualize the uncertainty
# in mu(i), the linear model of the mean. But actual predictions of heights depend also upon the
# stochastic definition in the first line. The Gaussian distribution on the first line tells us that
# the model expects observed heights to be distributed around mu, not right on top of it. And
# the spread around mu is governed by sigma. All of this suggests we need to incorporate sigma in the
# predictions somehow.

# Imagine simulating heights. For any unique weight value, you sample
# from a Gaussian distribution with the correct mean mu for that weight, using the correct
# value of sigma sampled from the same posterior distribution. If you do this for every sample
# from the posterior, for every weight value of interest, you end up with a collection of simulated
# heights that embody the uncertainty in the posterior as well as the uncertainty in the Gaussian likelihood.
sim.height <- sim( m4.3 , data=list(weight=weight.seq), n=1e5 )
str(sim.height)

# This matrix is much like the earlier one, mu, but it contains simulated heights, not distributions
# of plausible average height, mu.

# We can summarize these simulated heights in the same way we summarized the distributions of mu:
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )
height.mean <- apply (sim.height, 2, mean)
height.HPDI <- apply ( sim.height , 2 , HPDI , prob=0.89 )
#Now height.PI contains the 89% posterior prediction interval of observable 
#(according tothe model) heights, across the values of weight in weight.seq.

# Let’s plot everything we’ve built up: (1) the MAP line, (2) the shaded region of 89%
# plausible mu, and (3) the boundaries of the simulated heights the model expects.

# plot raw data
plot( height ~ weight , d2 , col=col.alpha(rangi2,0.5) )
# draw MAP line
lines( weight.seq , mu.mean )
# lines(weight.seq, height.mean) - just checking it's the same
# draw HPDI region for line
shade( mu.HPDI , weight.seq )
# draw PI region for simulated heights
shade( height.PI , weight.seq )

#The wide shaded region in the figure represents the
#area within which the model expects to find 89% of actual heights in the population, at each
#weight.

#4.5. Polynomial regression (with one predictor variable)
#In this context, “polynomial” means equations for mu(i) that add additional terms with squares, cubes, 
#and even higher powers of the predictor variable. There’s still only one predictor variable in the model, 
#so this is still a bivariate regression. But the definition of mu(i) has more parameters now

data(Howell1)
d <- Howell1

plot(d$height ~ d$weight)

# Fitting (polynomial) model to data

#1. Standardize
#Centre the variable (extract the mean) and divide by std
d$weight.s <- ( d$weight - mean(d$weight) )/sd(d$weight)

plot(d$height ~ d$weight.s)

#2. Fit

d$weight.s2 <- d$weight.s^2
m4.5 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b1*weight.s + b2*weight.s2 ,
    a ~ dnorm( 178 , 100 ) ,
    b1 ~ dnorm( 0 , 10 ) ,
    b2 ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d )

#3. Plot
#3.1 Calculate the mean relationship and the 89% intervals of the mean and the predictions
weight.seq <- seq( from=-2.2 , to=2 , length.out=30 )
pred_dat <- list( weight.s=weight.seq , weight.s2=weight.seq^2 )
mu <- link( m4.5 , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.89 )
sim.height <- sim( m4.5 , data=pred_dat )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

#3.2 Plot
plot( height ~ weight.s , d , col=col.alpha(rangi2,0.5) )
lines( weight.seq , mu.mean )
shade( mu.PI , weight.seq )
shade( height.PI , weight.seq )

#4.Cubic regression of weight
d$weight.s3 <- d$weight.s^3
m4.6 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b1*weight.s + b2*weight.s2 + b3*weight.s3 ,
    a ~ dnorm( 178 , 100 ) ,
    b1 ~ dnorm( 0 , 10 ) ,
    b2 ~ dnorm( 0 , 10 ) ,
    b3 ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d )

weight.seq <- seq( from=-2.2 , to=2 , length.out=30 )
pred_dat <- list( weight.s=weight.seq , weight.s2=weight.seq^2, weight.s3=weight.seq^3 )
mu <- link( m4.6 , data=pred_dat )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.89 )
sim.height <- sim( m4.6 , data=pred_dat )
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )


plot( height ~ weight.s , d , col=col.alpha(rangi2,0.5) )
lines( weight.seq , mu.mean )
shade( mu.PI , weight.seq )
shade( height.PI , weight.seq )

#Converting back to natural scale.
plot( height ~ weight.s , d , col=col.alpha(rangi2,0.5) , xaxt="n" )

at <- c(-2,-1,0,1,2)
labels <- at*sd(d$weight) + mean(d$weight)
axis( side=1 , at=at , labels=round(labels,1) )


#4.7 Practice

#4M1. Simulate observed heights from the prior (not the posterior) - model in book
sample_mu <- rnorm( 1e4 , 0 , 10 )
sample_sigma <- runif( 1e4 , 0 , 10 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma )
dens( prior_h )


#4M2. Translate the model into a map formula
# Won't actually work in practice, as there is no data for it
mapFormula <- alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(0,10),
    sigma ~ dunif(0,10)
      )

#4M4. A sample of students is measured for height each year for 3 years. After the third year, you want
# to fit a linear regression predicting height using year as a predictor. Write down the mathematical
# model definition for this regression, using any variable names and priors you choose. Be prepared to
# defend your choice of priors.



mp4.8 <- map(
  alist(
    height ~ dnorm(a+b*)
  )
)

#4H1.

data(Howell1) 
d <- Howell1

d$weight.c <- d$weight - mean(d$weight)

# Fit model
pModel1 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*weight ,
    a ~ dnorm( 139 , 100 ) ,
    b ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d )

precis( pModel1 , corr=TRUE )

plot( height ~ weight, data=d )
abline( a=coef(pModel1)["a"] , b=coef(pModel1)["b"])

# define sequence of weights to compute predictions for
# these values will be on the horizontal axis
weight.seq <- c(46.95, 43.72, 64.78, 32.59, 54.63)

# use link to compute mu for each sample from posterior
# and for each weight in weight.seq

muLink <- link(pModel1, data=data.frame(weight=weight.seq))

# or use the long form (of link)
post <- extract.samples(pModel1, n = 1e3) 
mu.link <- function(weight) post$a + post$b*log(weight)
#weight.seq <- (weight.seq - mean(weight.seq))/sd(weight.seq)
mu <- sapply( weight.seq , mu.link )


mu.mean <- apply( mu , 2 , mean )
mu.HPDI <- apply( mu , 2 , HPDI , prob=0.89 )

muLink.mean <- apply( muLink , 2 , mean )
muLink.HPDI <- apply( muLink , 2 , HPDI , prob=0.89 )

# Plot for sanity check
plot( height ~ weight , data=d , col=col.alpha(rangi2,0.5) )
# plot the MAP line, aka the mean mu for each weight
lines( weight.seq , muLink.mean )
# plot a shaded region for 89% HPDI
shade( muLink.HPDI , weight.seq )

# This only tells you the variation around the mean, says nothing about individual heights
# For that, you need to incorporate sigma in the predictions.

sim.height <- sim( pModel1 , data=list(weight=weight.seq) )
str(sim.height)

# Or use the long for (of sim)

post <- extract.samples(pModel1)
sim.height <- sapply( weight.seq , function(weight)
  rnorm(
    n=nrow(post) ,
    mean=post$a + post$b*weight ,
    sd=post$sigma ) )

height.HPDI <- apply( sim.height , 2 , HPDI , prob=0.89 )
height.mean <- apply( sim.height , 2 , mean )


# Repeat with new set of data (adults only)

# Filter data by age (i.e. adults only)
d2 <- d[d$age >= 18,]

d2$weight.c <- (d2$weight - mean(d2$weight))

# Replace d2$weight with d2$wight.c in the model fit:
pModel2 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*weight ,
    a ~ dnorm( 154 , 100 ) ,
    b ~ dnorm( 0 , 10 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d2)

precis( pModel2 , corr=TRUE )

plot( height ~ weight, data=d2 )
abline( a=coef(pModel2)["a"] , b=coef(pModel2)["b"])

sim.height2 <- sim( pModel2 , data=list(weight=weight.seq) )
str(sim.height)

height2.HPDI <- apply( sim.height2 , 2 , HPDI , prob=0.89 )
height2.mean <- apply( sim.height2 , 2 , mean )

# The second model appears to be a better fit (tighter HPDI).
# This on its own wouldn't be a good reason to prefer it over the first one ( - possibly overfitting)

# However, given the fact that the weights based on which the ind. heights need to be
# predicted all appear to be from adults (based on the data), it does make sense to 
# use this second model for our model-based predictions.



#4H2
data("Howell1")
d<-Howell1
rm(Howell1)
d2<- d[d$age < 18,]

#1.Fit a linear regression to these data, using map. Present and interpret the estimates. For every
# 10 units of increase in weight, how much taller does the model predict a child gets?
pModel3 <- map(
  alist(
    height ~ dnorm(a + b * weight, sigma),
    a ~ dnorm(100, 50),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ), 
  data = d2
)

precis(pModel3, corr = T)

# b = 2.72 => For every 10 units increase in weight, the model predicts a child will be 27.2 cm taller

#2. Plot the raw data, with height on the vertical axis and weight on the horizontal axis. Superimpose
# the MAP regression line and 89% HPDI for the mean. Also superimpose the 89% HPDI for
# predicted heights.

# Plot raw data with height on vertical axis and weight on horizontal axis
plot( height ~ weight, d2, col=col.alpha(rangi2,0.5) )

# Superimpose MAP regression line
abline( a=coef(pModel3)["a"] , b=coef(pModel3)["b"])

# Superimpose 89% HPDI for the mean
#mu <- link(pModel3, data=data.frame(weight=weight.seq))
# link doesn't work here. Why?!

weight.seq <- seq(round(min(d2$weight)), round(max(d2$weight)), by = 1)
post <- extract.samples(pModel3, n = 1e3) 
mu.link <- function(weight) post$a + post$b*weight
mu <- sapply( weight.seq , mu.link )

mu.mean <- apply( mu , 2 , mean )
mu.HPDI <- apply( mu , 2 , HPDI , prob=0.89 )

# plot the MAP line, aka the mean mu for each weight
lines(mu.mean, weight.seq)
# plot a shaded region for 89% HPDI
shade( mu.HPDI , weight.seq )

# Also superimpose the 89% HPDI for predicted heights.
sim.height <- sim(pModel3, data=data.frame(weight=weight.seq))
str(sim.height)

height.HPDI <- apply(sim.height, 2, HPDI, prob = .89)

shade( height.HPDI, weight.seq )

#3. What aspects of the model fit concern you? Describe the kinds of assumptions you would
# change, if any, to improve the model. You don’t have to write any new code. Just explain what the
# model appears to be doing a bad job of, and what you hypothesize would be a better model.

# The model's assumption that height and weight have the same relationship at every level.
# Specifically, this data suggest that there is a steep increase in height when weight is below ~30,
# after which height assimptotes. This likely corresponds to an age-related growth spur. Thus, 
# the includion of age as a predictor variable would likely offer a better model for the observed data.

#4H3.

#1. 
data("Howell1")
d <- Howell1
rm(Howell1)

pModel4 <- map(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*log2(weight) ,
    a ~ dnorm( 154 , 100 ) ,
    b ~ dnorm( 0 , 100 ) ,
    sigma ~ dunif( 0 , 50 )
  ) ,
  data=d)

precis(pModel4, corr = T)

plot( height ~ weight , d , col=col.alpha(rangi2,0.4) )


#2. 

# Use samples from the model to describe the model


# define sequence of weights to compute predictions for
# these values will be on the horizontal axis
weight.seq <- seq(1, 70, length.out = 150)
#weight.seq <- seq(round(min(log2(d$weight))), round(max(log2(d$weight))), by = 0.1)

# compute mu for each sample from posterior
# and for each weight in weight.seq
post <- extract.samples(pModel4, n = 1e3) 
mu.link <- function(weight) post$a + post$b*log2(weight)
#weight.seq <- (weight.seq - mean(weight.seq))/sd(weight.seq)
mu <- sapply( weight.seq , mu.link )

# Calculate mean and HPDI around mean based on samples from the model
mu.mean <- apply( mu , 2 , mean )
mu.HPDI <- apply( mu , 2 , HPDI , prob=0.97 )

# This only tells you the variation around the mean, says nothing about individual heights
# For that, you need to incorporate sigma in the predictions.
post <- extract.samples(pModel4, 1e3)
sim.height <- sapply( weight.seq , function(weight)
  rnorm(
    n=nrow(post) ,
    mean=post$a + post$b*log2(weight) ,
    sd=post$sigma ) )

# Calculate mean and HPDI around mean based on simulations based on the posterior
height.mean <- apply( sim.height , 2 , mean )
height.HPDI <- apply( sim.height , 2 , HPDI , prob=0.97 )


#Plot
# the predicted mean height as a function of weight (MAP line - mu for each weight in weight.seq)
lines( weight.seq , mu.mean )
# double check that the means from the simulations are the same?
lines( weight.seq , height.mean )

# plot a shaded region for 89% HPDI
shade( mu.HPDI , weight.seq )

# the 97% HPDI for the mean
shade( mu.HPDI , weight.seq )

# the 97% HPDI for predicted heights.
shade (height.HPDI ,weight.seq)







## 4H3 (from cavaunpeu)
data(Howell1)
d <- Howell1
trials <- 1e5

model <- map(
  alist(
    height ~ dnorm(mean = mu, sd = sigma),
    mu <- alpha + beta*log(weight),
    alpha ~ dnorm(mean = 178, sd = 100),
    beta ~ dnorm(mean = 0, sd = 100),
    sigma ~ dunif(min = 0, max = 50)
  ),
  data = d
)

# simulate mu then compute mean and hpdi
weight.seq <- seq(from = 1, to = 70, length.out = 100)
posterior.samples <- data.frame( mvrnorm(n = trials, mu = coef(model), Sigma = vcov(model)) )
mu.link <- function(weight) posterior.samples$alpha + posterior.samples$beta * log(weight)
mu <- sapply(X = weight.seq, FUN = mu.link)
mu.mean <- apply(X = mu, MARGIN = 2, FUN = mean)
mu.hpdi <- apply(X = mu, MARGIN = 2, FUN = HPDI, prob = .89)

# simulate heights then compute hpdi
height.link <- function(weight) rnorm(n = nrow(posterior.samples), mean = mu.link(weight), sd = posterior.samples$sigma)
height.samples <- sapply(X = weight.seq, FUN = height.link)
height.hpdi <- apply(X = height.samples, MARGIN = 2, FUN = HPDI, prob = .89)
height.mean <- apply(height.samples, 2, mean)

# plot results
plot(height ~ weight, data = d, col = col.alpha(rangi2, .4))
lines(x = weight.seq, y = mu.mean)
shade(object = mu.hpdi, lim = weight.seq)
shade(object = height.hpdi, lim = weight.seq)








