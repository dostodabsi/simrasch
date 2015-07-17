library('eRm')
library('lattice')
library('doParallel')


#' wrapper around RM function from eRm
#'
#' tryCatch until the simulated data matrix is neither ill-conditioned
#' nor has a participant with all 0 or all 1
#' @param n, items: numeric
#' @param model_sim: function to simulate Rasch data
#' @param ...: additional arguments to model_sim
#' @return returns the estimated eRm object
#' @examples
#' estimate(400, 30, sim.2pl, .5)
estimate <- function(n, items, model_sim, ...) {
  tryCatch({
    RM(model_sim(n, items, ...))
  }, warning = function(e) {
    estimate(n, items, model_sim, ...)
  }, error = function(e) {
    estimate(n, items, model_sim, ...)
  })
}

#' Get the itemfit statistics (infit, outfit, p-values)
#'
#' @param n, items, times: numeric
#' @param model_sim: function to simulate Rasch data
#' @param ...: additional arguments to model_sim
#' @return returns an object of class 'sim_res'
#' which is an array that stores #times #items x 3 matrices
#' @examples
#' get_itemfit(400, 30, sim.2pl, .4)
get_itemfit <- function(n, items, times, parallel, model_sim, ...) {

  if (parallel) {
    sim_res <- foreach(i = 1:times) %dopar% {
      rasch <- estimate(n, items, model_sim, ...)
      fits <- itemfit(person.parameter(rasch))
      cbind(fits$i.infitMSQ, fits$i.outfitMSQ,
            round(1 - pchisq(fits$i.fit, fits$i.df - 2), 3))
    }
    sim_res <- array(unlist(sim_res), dim = c(items, 3, times))
  } else {
    sim_res <- array(NA, dim = c(items, 3, times))
    for (i in 1:times) {
      rasch <- estimate(n, items, model_sim, ...)
      fits <- itemfit(person.parameter(rasch))
      sim_res[, , i] <- cbind(fits$i.infitMSQ, fits$i.outfitMSQ,
                              round(1 - pchisq(fits$i.fit, fits$i.df - 1), 3))
    }
  }

  dimnames(sim_res) <- list(1:items, c('in', 'out', 'p.val'), 1:times)
  attr(sim_res, 'n') <- n
  attr(sim_res, 'nitems') <- items
  attr(sim_res, 'class') <- 'sim_res'
  sim_res
}


#' Compute rejection of Rasch model based on Likelihood Ratio test after
#' items have been removed on the basis of item infit statistics
#'
#' @param n, items, times: numeric (single value)
#' @param model_sim: function to simulate Rasch data
#' @param ...: additional arguments to model_sim
#' @param cores: number of cores to run in parallel
#' @param cutoff: item infit statistic cutoff
#' @return returns an object of class 'sim_res'
#' which is an array that stores #times #items x 3 matrices
#' @examples
#' simrem(100, 20, sim.2pl, .4)
simrem <- function(n, items, times, model_sim, ...,
                   cores = 3, cutoff = c(.8, 1.2)) {

  cat('running in parallel\n') # always >:)
  registerDoParallel(cores = cores)

  low <- cutoff[1]
  high <- cutoff[2]
  pvals <- foreach(k = 1:times) %dopar% {
    fitted <- estimate(n, items, model_sim, ...)
    pp <- person.parameter(fitted)
    infit <- unname(itemfit(pp)$i.infitMSQ)
    remove <- which(infit < low | infit > high)
    if (length(remove) != 0) {
      reducedX <- fitted$X[, -remove]
      refitted <- RM(reducedX)
    } else {
      refitted <- fitted
    }
    LRtest(refitted)$pvalue
  }
  unlist(pvals)
}


#' Runs the main simulation
#'
#' @param times: numeric
#' @param n, items: numeric vectors
#' @param model_sim: function to simulate Rasch data
#' @param Sigma: factor correlation matrix for multidimensional data
#' @param weights: item loadings of the factor structure
#' @param parallel: if simulation should be run in parallel
#' @param cores: number of cores to run in parallel
#' @param ...: additional arguments to model_sim
#' @return returns a list that stores objects of class 'sim_res'
#' each array stores #times simulations of each possible n x items combination
#' @examples
#' sim(c(100, 400), c(10, 50), 1000, sim.2pl, .5, parallel = TRUE, cores = 3)
sim <- function(n, items, times, model_sim, ..., Sigma = NULL,
                weights = NULL, parallel = TRUE, cores = detectCores()) {
  if (parallel) {
    cat('running in parallel\n')
    registerDoParallel(cores = cores) # for parallel execution
  }

  g <- 1
  nlength <- length(n)
  nitems <- length(items)
  res <- rep(list(NULL), nlength * nitems)

  for (i in 1:nlength) {
    for (j in 1:nitems) {

      if (!is.null(Sigma) && !is.null(weights)) {
        factors <- ncol(Sigma)
        wmat <- apply(replicate(items[j] / factors, weights), 2, rbind)
        res[[g]] <- get_itemfit(n[i], items[j], times,
                      parallel = parallel, sim.xdim, Sigma, wmat)
      } else {
        res[[g]] <- get_itemfit(n[i], items[j], times,
                      parallel = parallel, model_sim, ...)
      }

      cat('done\n')
      g <- g + 1
    }
  }
  attr(res, 'simtype') <- c(as.character(substitute(model_sim)),
                            as.character(substitute(...)))
  res
}


#' Compute the items suggested to be removed / revised based on item infit
#'
#' @param sim_res: simulation results
#' @param cutoff: itemfit cutoff
#' @return returns a list with alpha / beta errors across all
#' person x item combinations aggregated over all matrices, and the
#' number of items removed (for each matrix)
#' @examples
#' compute_errors(sim_res)
compute_errors <- function(sim_res, cutoff = c(.8, 1.2)) {
  low <- cutoff[1]
  high <- cutoff[2]
  error_res <- list()
  nr_sims <- length(sim_res)

  for (i in 1:nr_sims) {
    pvals <- 0
    cursim <- sim_res[[i]]

    times <- dim(cursim)[3]
    nr_items <- attr(cursim, 'nitems')
    normalizer <- times * nr_items
    rm_items <- rep(NA, times)

    for (j in 1:times) {
      rm_items_mat <- 0
      for (k in 1:nr_items) {
        instat <- cursim[k, , j][1]
        pval <- cursim[k, , j][3]

        if (pval <= .05) {
          pvals <- pvals + 1
        }

        if (instat < low || instat > high) {
          rm_items_mat <- rm_items_mat + 1
        }
      }
      rm_items[j] <- rm_items_mat
    }
    error_res[[i]] <- list('rm_percent' = sum(rm_items) / normalizer,
                           'n' = attr(cursim, 'n'), 'nr_items' = nr_items,
                           'rm_items' = rm_items, 'pval' = pvals / normalizer)
  }
  error_res
}


#' Store the results of the compute_err function in a data.frame
#'
#' @param sim_res: result of the simulation (list of arrays)
#' @param cutoff: cutoff used for infit statistics (commonly .8 - 1.2)
#' @return returns a data.frame with four columns: the number of participants (n),
#' the number of items (items), and the percent (!) that pvalues and infit
#' statistics respectively would reject an item. Thus note that the column
#' p-value and infit statistics are **not** p-values nor infit statistics.
#' It is a percent which is based on the number of rejections of items due to
#' standardized p-values  < .05 or infit statistics outside a certain cutoff
#' @examples
#' summarize(err_res)
summarize <- function(sim_res, cutoff = c(.8, 1.2)) {
  err_res <- compute_errors(sim_res, cutoff = cutoff)
  pick <- function(sim) list('sims' = length(sim$rm_items),
                             'n' = sim$n, 'items' = sim$nr_items,
                             'pval' = round(sim$pval, 3),
                             'infit' = round(sim$rm_percent, 3),
                             'simtype' = attr(sim_res, 'simtype')[1],
                             'violation' = attr(sim_res, 'simtype')[2])
  df <- data.frame(t(sapply(err_res, pick)))
  df <- data.frame(apply(df, 2, unlist))
  df$infit <- as.numeric(as.character(df$infit))
  df$items <- factor(df$items, sort(as.numeric(levels(df$items))))
  if (all(is.na(df$violation))) {
    df$violation <- ifelse(df$simtype == 'sim.rasch', 'rasch', 'mult')
  }
  df
}


#' Visualizes the \% items removed / Rasch rejected given a simulation result
#'
#' @param res_sum: summarized simulation result or processed LR simulation
#' @param ylab: label of the y-axis
#' @param ...: additional arguments to xyplot
#' @return plots a trellis graph visualising the results
#' @examples
#' visualize(res_sum)
visualize <- function(res_sum, ylab = '% removed items', ...) {
  if (all(is.na(res_sum$violation))) res_sum$violation <- 'none'
  if ('infit' %in% names(res_sum)) res_sum$reject <- res_sum$infit

  xyplot(reject ~ items | n, groups = violation, data = res_sum,
         type = c('p', 'g'),
         xlab = list(label = 'item size', cex = 2),
         ylab = list(label = ylab, cex = 2),
         auto.key = list(columns = 2, cex = 2),
         par.settings = simpleTheme(pch = 21, cex = 2),
         par.strip.text = list(cex = 1.8),
         axis.text = list(cex = 2),
         scales = list(cex = 1.6),
         panel = function(...) {
           panel.abline(h = .05, lty = 'dotted', col = 'black', lwd = 2)
           panel.xyplot(...)
       }, ...)
}
