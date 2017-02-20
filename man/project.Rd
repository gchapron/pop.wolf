\name{project}
\alias{project}
\alias{C_montecarlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wolf population projections
}
\description{
Run stochastic wolf population projections.
}
\usage{
project(years, runs,
     surv_pup, surv_sub, surv_vag, surv_adt,
     dispers_weib_shape, dispers_weib_scale,
     settl_weib_shape, settl_weib_scale,
     pair1breed, litter_size,
     pop_initial, pop_quota, seed)
}
%- maybe also 'usage' for other objects documented here.

\arguments{

  \item{years}{
A number: number of years to project the population.
}

  \item{runs}{
A number: number of times (or Monte Carlo runs) to project the population.
}

  \item{surv_pup}{
A vector: average and standard deviation of annual survival of pups from 0 to 5 months.
}

  \item{surv_sub}{
A vector: average and standard deviation of annual survival of subadults (>5 month-old non breeding animals in packs).
}

  \item{surv_vag}{
A vector: average and standard deviation of annual survival of vagrants.
}

  \item{surv_adt}{
A vector: average and standard deviation of survival of resident adults.
}

  \item{dispers_weib_shape}{
A vector: average and standard deviation of the shape of a Weibul distributed variable describing the time in months until dispersal.
}

  \item{dispers_weib_scale}{
A vector: average and standard deviation of the scale of a Weibul distributed variable describing the time in months until dispersal.
}

  \item{settl_weib_shape}{
A vector: average and standard deviation of the shape of a Weibul distributed variable describing the time in months from dispersal until settlement.
}

  \item{settl_weib_scale}{
A vector: average and standard deviation of the scale of a Weibul distributed variable describing the time in months from dispersal until settlement.
}

  \item{pair1breed}{
Number of times (or Monte Carlo runs) to project the population.
}

  \item{litter_size}{
A vector: average and standard deviation of number of pups in litter.
}

  \item{pop_initial}{
(optional) A list: with elements packs and vagrants. \code{packs} is an array representing the composition of the initial population, rows are packs, columns indicate the age of the breeding male, the age of the breeding female, and the pack size. \code{vagrants} is a number representing the number of vagrants in the initial population.
}

  \item{pop_quota}{
(optional) An array: number of animals killed every month (rows) from the following categories (columns): pairs, alphas, vagrants, subadults, unknown. Removing a pair kills 2 animals but is coded as 1.
}

  \item{seed}{
(optional) A number: seed of the random number generator.
}


}
\details{
Run stochastic wolf population projections with an Individual-Based Model (IBM) compiled in C.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{runs}{a 3-dimensional array of numbers of individuals with dimension c(years, statistics, runs)}

  \item{individuals}{a 2-dimensional array of individuals events}

  \item{parameters}{a list of parameters of the projection}
%% ...
}

\examples{

# Number of years to project
years <- 15

# Initial population
pop_initial <- list()
# 10 packs of size 5 with alpha male aged 2 and alpha female aged 3
pop_initial$packs <- matrix(rep(c(2,3,5),10), ncol=3, nrow=10, byrow=TRUE)
# 5 vagrants
pop_initial$vagrants <- 5

# Hunting quota
pop_quota <- matrix(0, nrow=12*years+1, ncol=5)
# Hunt 4 alphas every March
pop_quota[1+seq(3,years*12,12),2] <- 4

projection <- project(
	years = years,
	runs = 100,
	surv_pup = c(0.65, 0.1),
	surv_sub = c(0.74, 0.13),
	surv_vag = c(0.42, 0.1),
	surv_adt = c(0.82, 0.03),
	dispers_weib_shape = c(1.419, 0.154),
	dispers_weib_scale = c(4.687, 0.528),
	settl_weib_shape = c(1.353, 0.210),
	settl_weib_scale = c(4.606, 0.772),
	pair1breed = c(0.79, 0.05),
	litter_size = c(5, 0.52),
	pop_initial = pop_initial,
	pop_quota = pop_quota,
	seed = 1
)
}

