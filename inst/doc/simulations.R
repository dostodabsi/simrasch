## ------------------------------------------------------------------------
library('eRm')
library('lattice')
library('simrasch')

items <- c(12, 24, 60, 120) # number of items
n <- c(150, 300, 500, 1000) # number of participants
weights <- diag(.5, 4, 4) # item loadings on each factor
sigma <- matrix(c(1, 0.6, 0.5, 0.4, # factor loadings
                  0.6, 1, 0.7, 0.4, # for a theoretically plausible
                  0.5, 0.7, 1, 0.3, # 4-dimensional model of scientific reasoning
                  0.4, 0.4, 0.3, 1), 4)

## ---- eval = FALSE-------------------------------------------------------
#  # no violation
#  sim_rasch <- sim(n, items, 1000, sim.rasch, parallel = TRUE, cores = 4)
#  
#  # violation of equal item discrimination
#  sim_2pl1 <- sim(n, items, 1000, sim.2pl, .1, parallel = TRUE, cores = 4)
#  sim_2pl3 <- sim(n, items, 1000, sim.2pl, .3, parallel = TRUE, cores = 4)
#  sim_2pl5 <- sim(n, items, 1000, sim.2pl, .5, parallel = TRUE, cores = 4)
#  
#  # assing fit of Rasch model using Andersen's LR test
#  # after removing items based on infit (example call for items = 120, n = 1000)
#  sim_rem <- simrem(1000, 120, 1000, sim.2pl, .5, cores = 4)
#  
#  # violation of unidimensionality
#  sim_mult <- sim(n, items, 1000, sim.xdim,
#                  Sigma = sigma, weights = weights, parallel = TRUE, cores = 4)

## ---- message = FALSE----------------------------------------------------
data(sim_2pl1, sim_2pl3, sim_2pl5, sim_rasch, sim_mult, sim_removal)

## ------------------------------------------------------------------------
sm_rasch <- summarize(sim_rasch)
sm_rasch

## ---- fig.width = 7, fig.height = 7, warning = FALSE---------------------
processed <- data.frame(violation = rep(c('pval', 'infit'), each = nrow(sm_rasch)),
                        reject = c(sm_rasch$pval, sm_rasch$infit), n = rep(sm_rasch$n, 2),
                        items = rep(sm_rasch$items, 2))
visualize(processed, ylim = c(0, .3))

## ------------------------------------------------------------------------
dat <- sim.2pl(1000, 60, discrim = .5)
mrasch <- RM(dat)
pp <- person.parameter(mrasch)
eRm::itemfit(pp)

## ---- fig.width = 7, fig.height = 7--------------------------------------
sm_pl <- rbind(summarize(sim_2pl1), summarize(sim_2pl3), summarize(sim_2pl5))
visualize(sm_pl, ylim  = c(0, .3))

## ---- message = FALSE----------------------------------------------------
library('dplyr')
processed <- sim_removal %>%
             group_by(n, items, model, violation) %>%
             summarise(reject = mean(pvals <= .05)) %>% data.frame

processed

## ---- fig.width = 7, fig.height = 7--------------------------------------
visualize(processed, ylab = '% rejected Rasch')

## ---- message = FALSE----------------------------------------------------
library('mirt')

mrasch <- mirt(dat, 1, itemtype = 'Rasch', verbose = FALSE)
mbirn <- mirt(dat, 1, itemtype = '2PL', verbose = FALSE)

anova(mrasch, mbirn)

## ------------------------------------------------------------------------
items <- 60
wmat <- apply(replicate(items / 4, weights), 2, rbind) # four factors
dat <- sim.xdim(1000, items, Sigma = sigma, weightmat = wmat)

mrasch <- RM(dat)
pp <- person.parameter(mrasch)
eRm::itemfit(pp)

## ------------------------------------------------------------------------
sm_mult <- simrasch::summarize(sim_mult)

## ---- warning = FALSE, fig.width = 7, fig.height = 7---------------------
visualize(sm_mult, ylim = c(0, .3))

## ------------------------------------------------------------------------
MLoef(mrasch, rep(1:4, times = items / 4))

## ---- eval = FALSE-------------------------------------------------------
#  # takes over one hour to run
#  res <- NPtest(dat, n = 500, method = 'MLoef', splitcr = rep(1:4, items / 4))

## ------------------------------------------------------------------------
readRDS('NP_res.RDS') # load the results from disk

## ------------------------------------------------------------------------
model <- 'f1=1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57
          f2=2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58
          f3=3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59
          f4=4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60
          cov=f1*f2*f3*f4'

fitm <- mirt(dat, 1, itemtype = 'Rasch', verbose = FALSE)
#fit2pl <- mirt(dat, 4, itemtype = '2PL')

## ------------------------------------------------------------------------
fscores(fitm, returnER = TRUE)

