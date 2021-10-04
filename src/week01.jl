using Distributions
"""
    function week01_1()
---
# Description
Suppose the globe tossing data (Chapter 2) had turned out to be 4 water
in 15 tosses.  Construct the posterior distribution, using grid approximation.
Use the same flat prior as in the book.

# My understanding
Suppose the true parameter `p` is uniformly distributed.
Then for "gridded" `p`, the probability fo the observations are different.

# R codes
```R
## Problem 1

# define grid
p_grid <- seq( from=0 , to=1 , length.out=1000 )
# define prior
prior <- rep( 1 , 1000 )
# compute likelihood at each value in grid
likelihood <- dbinom( 4 , size=15 , prob=p_grid )
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="l" , xlab="proportion water" , ylab="posterior probability" )
```
"""
function week01_1()
    grid = 0:.001:1
    ll = pdf.(Binomial.(15, grid), 4) # likelihood
    post = ll ./ sum(ll)
    plot(grid, post, xlabel = "proportion water", ylabel = "posterior probability", legend = false)
end

"""
    function week01_2()
---
# Description
Start over in 1, but now use a prior that is zero below p = 0.5 and a con-
stant above p = 0.5. This corresponds to prior information that a majority
of the Earth’s surface is water. What difference does the better prior make?

# R codes
```R
# define grid
p_grid <- seq( from=0 , to=1 , length.out=1000 )
# define prior
prior <- c( rep( 0 , 500 ) , rep( 2 , 500 ) )
# compute likelihood at each value in grid
likelihood <- dbinom( 4 , size=15 , prob=p_grid )
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="l" , xlab="proportion water" , ylab="posterior probability" )
```

#implementation
We can just define the grid of `0.5:.001:1`.  All other likelihood/posteriors are zero.
"""
function week01_2()
    grid = 0.5:.001:1
    ll = pdf.(Binomial.(15, grid), 4)
    posterior = ll ./ sum(ll)
    plot(grid, posterior,
         xlabel = "proportion water",
         ylabel = "posterior probability",
         legend = false)
end

"""
    function week01_3()
---
# Description
For the posterior distribution from 2, compute 89% percentile and HPDI
intervals. Compare the widths of these intervals. Which is wider? Why? If
you had only the information in the interval, what might you misunderstand
about the shape of the posterior distribution?

# R codes
```R
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )
plot( samples , ylim=c(0,1) , xlab="samples" , ylab="proportion water" )
PI( samples , 0.89 )
HPDI( samples , 0.89)
```
"""
function week01_3()
    @info "I skipped this assignment"
end

