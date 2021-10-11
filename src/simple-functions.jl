"""
    function PI(v::Vector{T}, α = 0.11) where {T <: Real}
---
# Description
Return the `P`ercentile `I`nterval of an array.

# Notes:
- `v`: is the posterior distribution and is corresponding 0 → step → 1.
- `α`: two-tail probability.

# R codes
```R
PCI <- function( samples , prob=0.89 ) {
    x <- sapply( prob , function(p) {
        a <- (1-p)/2
        quantile( samples , probs=c(a,1-a) )
    } )
    # now order inside-out in pairs
    n <- length(prob)
    result <- rep(0,n*2)
    for ( i in 1:n ) {
        low_idx <- n+1-i
        up_idx <- n+i
        # lower
        result[low_idx] <- x[1,i]
        # upper
        result[up_idx] <- x[2,i]
        # add names
        a <- (1-prob[i])/2
        names(result)[low_idx] <- concat(round(a*100,0),"%")
        names(result)[up_idx] <- concat(round((1-a)*100,0),"%")
    }
    return(result)
}
PI <- PCI
```
"""
function PI(b, p, α = .11)
    @info "Under construction"
end

"""
    function HPDI(b, p, α = .11)
---
# Description
`H`ighest `P`osterior `D`ensity `I`nterval (HDPI).
That is the narrorwest interval containing the specified probability mass.

The algorithm used in the course is to sample 1e4 samples `grids` with the
`posterior` as probability. Then sort the boundaries such that the probability
between the left and right boundary is `1 - α`.

My algorithm here is to find the peak of `posterior` first.
Then add left or right `posterior` probability, if its bigger.
Until the sum is greater than `1 - α`.

- `b`: the grid values, of which to find boundaries.
- `p`: the posterior probability of values in `b`.
- `α`: two tailed probability.

This algorithm can go wrong if there are more than one peaks.
It also need the length of `b` to be long enough.
The algorithm was not optimized, or considered for many pitfalls.

# R codes
```R
HPDI <- function( samples , prob=0.89 ) {
    # require(coda)
    coerce.list <- c( "numeric" , "matrix" , "data.frame" , "integer" , "array" )
    if ( inherits(samples, coerce.list) ) {
        # single chain for single variable
        samples <- coda::as.mcmc( samples )
    }
    x <- sapply( prob , function(p) coda::HPDinterval( samples , prob=p ) )
    # now order inside-out in pairs
    n <- length(prob)
    result <- rep(0,n*2)
    for ( i in 1:n ) {
        low_idx <- n+1-i
        up_idx <- n+i
        # lower
        result[low_idx] <- x[1,i]
        # upper
        result[up_idx] <- x[2,i]
        # add names
        names(result)[low_idx] <- concat("|",prob[i])
        names(result)[up_idx] <- concat(prob[i],"|")
    }
    return(result)
}
```
"""
function HPDI(b, p, α = .11)
    sum(p) < 1 - α && return b[1], b[end]
    s, c = findmax(p)
    l, r = c, c
    while s < 1 - α
        if l == 1
            r += 1
            s += p[r]
        elseif r == length(p)
            l -= 1
            s += p[l]
        elseif p[l] > p[r]
            l -= 1
            s += p[l]
        else
            r += 1
            s += p[r]
        end
    end
end

#=
Test codes of week 01
# Statistical Rethinking Winter 2020/2021
# Discussion Seminar 1
# Points
# (1) Procedural questions
# (2) Software setup problems
# (3) Homework review - solutions and broader concepts
# (4) Prepare for next week - Chapter 4!

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

## Problem 2

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

## Problem 3

samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )
plot( samples , ylim=c(0,1) , xlab="samples" , ylab="proportion water" )
PI( samples , 0.89 )
HPDI( samples , 0.89)
=#
