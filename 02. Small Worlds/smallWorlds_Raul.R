#### Autumn 2018 ########################################################################################

### Grid approximation (RCode 2.3)

# define grid
p_grid <- seq( from=0 , to=1 , length.out=20 )

# define prior
prior <- rep(1,20)
#prior <- ifelse( p_grid < 0.5 , 0 , 1 )
#prior <- exp( -5*abs( p_grid - 0.5 ) )

# compute LH at each value in grid
lh <- dbinom(6, size = 9, prob = p_grid)

# compute product of likelihood and prior
unstd.posterior <- lh * prior

# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

# plot 
plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )


### Quadratic approximation (RCode 2.6)

# To compute the quadratic approximation for the globe tossing data, we’ll use a tool in
# the rethinking package: map. The abbreviation MAP stands for maximum a posteriori,
# which is just a fancy Latin name for the mode of the posterior distribution.

# To use map, you provide a formula, a list of data, and a list of start values for the parameters.
# The formula defines the likelihood and prior.

library(rethinking)
globe.qa <- map(
  alist(
    w ~ dbinom(9,p) , # binomial likelihood
    p ~ dunif(0,1) # uniform prior
  ) ,
  data=list(w=6) )
# display summary of quadratic approximation
precis( globe.qa )

# You can read this kind of approximation like: Assuming the posterior is Gaussian, it is maximized
# at 0.67, and its standard deviation is 0.16.


# analytical calculation
w <- 24
n <- 36
curve( dbeta( x , w+1 , n-w+1 ) , from=0 , to=1 )
# quadratic approximation
curve( dnorm( x , 0.67 , 0.16 ) , lty=2 , add=TRUE )

# PRACTICE

#2E4. We can be 70% certain that the next toss will be water. 

#2M1.Recall the globe tossing model from the chapter. Compute and plot the grid approximate
# posterior distribution for each of the following sets of observations. In each case, assume a uniform
# prior for p.

# (1) W,W,W
p_grid <- seq( from=0 , to=1 , length.out=20 )
prior <- rep(1,20)
lh <- dbinom(3, size = 3, prob = p_grid)
unstd.posterior <- lh * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )

# (2) W,W,W,L
p_grid <- seq( from=0 , to=1 , length.out=20 )
prior <- rep(1,20)
lh <- dbinom(3, size = 4, prob = p_grid)
unstd.posterior <- lh * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )

# (3) L,W,W,L,W,W,W
p_grid <- seq( from=0 , to=1 , length.out=20 )
prior <- rep(1,20)
lh <- dbinom(5, size = 7, prob = p_grid)
unstd.posterior <- lh * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )

# 2M2. Now assume a prior for p that is equal to zero when p < 0:5 and is a positive constant when
# p >= 0:5. Again compute and plot the grid approximate posterior distribution for each of the sets of
# observations in the problem just above.

# (1) W,W,W
p_grid <- seq( from=0 , to=1 , length.out=20 )
prior <- rep(1,20)
prior <- ifelse(p_grid<0.5, 0, 1)
lh <- dbinom(3, size = 3, prob = p_grid)
unstd.posterior <- lh * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )

# (2) W,W,W,L
p_grid <- seq( from=0 , to=1 , length.out=20 )
prior <- rep(1,20)
prior <- ifelse(p_grid<0.5, 0, 1)
lh <- dbinom(3, size = 4, prob = p_grid)
unstd.posterior <- lh * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )

# (3) L,W,W,L,W,W,W
p_grid <- seq( from=0 , to=1 , length.out=20 )
prior <- rep(1,20)
prior <- ifelse(p_grid<0.5, 0, 1)
lh <- dbinom(5, size = 7, prob = p_grid)
unstd.posterior <- lh * prior
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )

#2M3. Suppose there are two globes, one for Earth and one for Mars. The Earth globe is 70% covered
# in water. The Mars globe is 100% land. Further suppose that one of these globes—you don’t know
# which—was tossed in the air and produced a “land” observation. Assume that each globe was equally
# likely to be tossed. Show that the posterior probability that the globe was the Earth, conditional on
# seeing “land” (Pr(Earthjland)), is 0.23.

prior <- c(.5, .5)
likelihood <- c(.3, 1)
unstandardized.posterior <- prior * likelihood
posterior <- unstandardized.posterior / sum(unstandardized.posterior)
round( posterior[1], 2) == .23
# See notebook for handwritten solution

#2M4. Suppose you have a deck with only three cards. Each card has two sides, and each side is either
# black or white. One card has two black sides. The second card has one black and one white side. The
# third card has two white sides. Now suppose all three cards are placed in a bag and shuffled. Someone
# reaches into the bag and pulls out a card and places it flat on a table. A black side is shown facing up,
# but you don’t know the color of the side facing down. Show that the probability that the other side is
# also black is 2/3. Use the counting method (Section 2 of the chapter) to approach this problem. This
# means counting up the ways that each card could produce the observed data (a black side facing up on the table).

ways <- c(2,1,0)
plau <- ways/sum(ways)
plau[1]



card.1.likelihood <- 2
card.2.likelihood <- 1
card.3.likelihood <- 0
likelihood <- c(card.1.likelihood, card.2.likelihood, card.3.likelihood)
prior <- rep(x = 1, length = length(likelihood))
unstandardized.posterior <- prior * likelihood
posterior <- unstandardized.posterior / sum(unstandardized.posterior)

#2M5. Now suppose there are four cards: B/B, B/W, W/W, and another B/B. Again suppose a card is
# drawn from the bag and a black side appears face up. Again calculate the probability that the other
# side is black.

ways <- c(2,1,0,2)
plau <- ways/sum(ways)

card.1.likelihood <- 2
card.2.likelihood <- 1
card.3.likelihood <- 0
card.4.likelihood <- 2
likelihood <- c(card.1.likelihood, card.2.likelihood, card.3.likelihood, card.4.likelihood)
prior <- rep(x = 1, length = length(likelihood))
unstandardized.posterior <- prior * likelihood
posterior <- unstandardized.posterior / sum(unstandardized.posterior)

# the probability the other size is black is equal to the probability that we've drawn card 1 or 4
posterior[1] + posterior[4]


#2M6. Imagine that black ink is heavy, and so cards with black sides are heavier than cards with white
# sides. As a result, it’s less likely that a card with black sides is pulled from the bag. So again assume
# there are three cards: B/B, B/W, and W/W. After experimenting a number of times, you conclude that
# for every way to pull the B/B card from the bag, there are 2 ways to pull the B/W card and 3 ways to
# pull the W/W card. Again suppose that a card is pulled and a black side appears face up. Show that
# the probability the other side is black is now 0.5. Use the counting method, as before.

ways <- c(2,2,0)
plau <- ways/sum(ways)

card.1.likelihood <- 2
card.2.likelihood <- 1
card.3.likelihood <- 0
likelihood <- c(card.1.likelihood, card.2.likelihood, card.3.likelihood)
prior <- c(1, 2, 3)
unstandardized.posterior <- prior * likelihood
posterior <- unstandardized.posterior / sum(unstandardized.posterior)

# the probability the other size is black is equal to the probability that we've drawn card 1
posterior[1] == .5

#2M7. Assume again the original card problem, with a single card showing a black side face up. Before
# looking at the other side, we draw another card from the bag and lay it face up on the table. The face
# that is shown on the new card is white. Show that the probability that the first card, the one showing
# a black side, has black on its other side is now 0.75. Use the counting method, if you can. Hint: Treat
# this like the sequence of globe tosses, counting all the ways to see each observation, for each possible
# first card.

# !!! FIGURE THIS OUT !!!
card.1.2.likelihood <- 2
card.2.1.likelihood <- 0
card.1.3.likelihood <- 4
card.3.1.likelihood <- 0
card.2.3.likelihood <- 2
card.3.2.likelihood <- 0

likelihood <- c(card.1.2.likelihood, card.2.1.likelihood, card.1.3.likelihood, card.3.1.likelihood, card.2.3.likelihood, card.3.2.likelihood)
prior <- rep(x = 1, length = length(likelihood))
unstandardized.posterior <- prior * likelihood
posterior <- unstandardized.posterior / sum(unstandardized.posterior)

# the probability that the other side of the first card is black is equal to the probability that the first card is card 1,
# which equals the probability that the sequence we've chosen is either (1, 2), or (1, 3)
posterior[1] + posterior[3] == .75

## THE REST IS IN THE NOTEBOOK!



#### WINTER 2019 ########################################################################################
#In NOTEBOOK









