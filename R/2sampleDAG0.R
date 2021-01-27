#' Implementation of the parallel-tempered MCMC for one-sample DAG0
#'
#' @param data0 a matrix containing data for the control group.
#' @param data1 a matrix containing data for the case/treatment group.
#' @param starting a list with each tag corresponding to a parameter name. 
#' Valid tags are 'E0', 'D', 'E1', 'U', 'zeta', 'eta' 'alpha0', 'alpha1', 'beta0', 'beta1', 'delta0', 'delta1', 'gamma0', 'gamma1', 'psi0', 'psi1', 'tau', and 'rho'. 
#' The value portion of each tag is the parameters' starting values for MCMC.
#' @param tuning a list with each tag corresponding to a parameter name. 
#' Valid tags are 'E', 'E_rev', 'U', 'zeta', 'eta', 'alpha0', 'alpha1', 'beta0', 'beta1', delta0', 'delta1', 'gamma0', 'gamma1', 'psi0', and 'psi1'. 
#' The value portion of each tag defines the standard deviations of Normal proposal distributions for Metropolis sampler for each parameter. 
#' The tag 'E' is for the birth or death update schemes for E = \{E_0, E_1\}, while 'E_rev' is for the reversal schemes.
#' The value portion of both tags should consist of the standard deviations of Gaussian proposals for alpha, beta, delta, gamma, zeta, and eta for an update of E.
#' Also, the value portion of the tag 'U' should include the standard deviations of Gaussian proposals for zeta and eta for an update of U.
#' @param priors a list with each tag corresponding to a parameter name. 
#' Valid tags are 'psi', 'tau', and 'rho'. The value portion of each tag defines the hyperparameters for two-sample DAG0.
#' @param n_sample the number of MCMC iterations.
#' @param n_burnin the number of burn-in samples.
#' @param n_chain the number of Markov chains for the parallel tempering.
#' @param prob_swap the probability of a swapping step.
#' @param temperature a sequence of increasing temperatures. Should have length of n_chain. 
#' Default value is T_l = (1.5)^\{(l - 1) / (n_chain - 1)\} for l = 1, ..., n_chain.
#' @param do_par if TRUE, the parallel computation is used to update n_chain Markov chains in parallel. 
#' The parallel computation is only available on Unix-alike platforms as we create the worker process by forking. 
#' @param n_core the number of cores to be used for the parallel computation when do_par = TRUE. 
#' @param verbose if TRUE, progress of the sampler is printed to the screen. Otherwise, nothing is printed to the screen.
#' @param n_report the interval to report MCMC progress.
#'
#' @return An object of class 1sampleDAG0, which is a list with the tags 'samples' and 'acceptance'. 
#' The value portion of the tag 'samples' gives MCMC samples from the posterior distribution of the parameters for two-sample DAG0.
#' The value portion of the tag 'acceptance' shows the Metropolis sampling acceptance percents for each parameter. 
#' @export
#'
#' @examples
#' To be written.
mcmc_2sampleDAG0 = function(data0, data1, starting, tuning, priors, n_sample = 5000, n_burnin = 2500, 
                            n_chain = 10, prob_swap = 0.1, temperature = NULL, do_par = FALSE, n_core = n_chain,
                            verbose = TRUE, n_report = 100)
{
   # sample size and dimension
   n0 = nrow(data0)
   n1 = nrow(data1)
   p  = ncol(data0)
   
   # default temperatures for the parallel tempering 
   if (is.null(temperature))
      temperature = exp(seq(0, log(1.5), length.out = n_chain))

   # sd's of the Metropolis sampler Normal proposal distribution for each chain
   list_tuning = list()
   for (m in 1 : n_chain)
   {
      list_tuning[[m]] = list()
      list_tuning[[m]]$E      = tuning$E
      list_tuning[[m]]$E_rev  = tuning$E_rev
      list_tuning[[m]]$U      = tuning$U
      list_tuning[[m]]$zeta   = tuning$zeta   * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$eta    = tuning$eta    * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$alpha0 = tuning$alpha0 * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$alpha1 = tuning$alpha1 * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$beta0  = tuning$beta0  * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$beta1  = tuning$beta1  * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$delta0 = tuning$delta0 * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$delta1 = tuning$delta1 * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$gamma0 = tuning$gamma0 * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$gamma1 = tuning$gamma1 * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$psi0   = tuning$psi0   * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$psi1   = tuning$psi1   * temperature[n_chain]^((m - 1) / (n_chain - 1))
   }
   
   # calculate logitPi and logMu for each group with starting values
   starting$logitPi0 = tcrossprod(data0, starting$alpha0) + matrix(starting$delta0, n0, p, byrow = TRUE)
   starting$logitPi1 = tcrossprod(data1, starting$alpha1) + matrix(starting$delta1, n1, p, byrow = TRUE)
   starting$logMu0   = tcrossprod(data0, starting$beta0) + matrix(starting$gamma0, n0, p, byrow = TRUE)
   starting$logMu1   = tcrossprod(data1, starting$beta1) + matrix(starting$gamma1, n1, p, byrow = TRUE)
   
   # initialize parameters
   param = list()
   for (m in 1 : n_chain)
   {
      param[[m]] = starting
   }
   
   # initialize MCMC samples for the cold chain
   mcmc_samples = list()
   mcmc_samples$E0     = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$D      = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$E1     = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$U      = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$zeta   = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$eta    = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$alpha0 = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$alpha1 = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$beta0  = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$beta1  = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$delta0 = matrix(NA, p, n_sample)
   mcmc_samples$delta1 = matrix(NA, p, n_sample)
   mcmc_samples$gamma0 = matrix(NA, p, n_sample)
   mcmc_samples$gamma1 = matrix(NA, p, n_sample)
   mcmc_samples$psi0   = matrix(NA, p, n_sample)
   mcmc_samples$psi1   = matrix(NA, p, n_sample)
   mcmc_samples$tau    = matrix(NA, 6, n_sample)
   mcmc_samples$rho    = matrix(NA, 3, n_sample)
   
   # initialize acceptance indicators for the cold chain
   acceptance = list()
   acceptance$zeta   = acceptance$eta    = array(NA, dim = c(p, p, n_sample))
   acceptance$alpha0 = acceptance$alpha1 = array(NA, dim = c(p, p, n_sample))
   acceptance$beta0  = acceptance$beta1  = array(NA, dim = c(p, p, n_sample))
   acceptance$delta0 = acceptance$delta1 = matrix(NA, p, n_sample)
   acceptance$gamma0 = acceptance$gamma1 = matrix(NA, p, n_sample)
   acceptance$psi0   = acceptance$psi1   = matrix(NA, p, n_sample)
   
   # run MCMC iterations
   if (do_par)
   {
      cluster = makeCluster(n_core, type = "FORK", setup_timeout = 0.5)
      registerDoParallel(cluster)
   }
   
   for (t in 1 : n_sample)
   {
      if (runif(1) > prob_swap)
      {
         # perform one-step update for all chains
         if (do_par)          # update chains in parallel
         {
            out = foreach(m = 1 : n_chain, .packages = "igraph") %dorng% {
               gibbs_2sampleDAG0(param[[m]], list_tuning[[m]], priors, data0, data1, temperature[m])
            }
         } else               # update chains sequentially
         {
            out = foreach(m = 1 : n_chain, .packages = "igraph") %do% {
               gibbs_2sampleDAG0(param[[m]], list_tuning[[m]], priors, data0, data1, temperature[m])
            }
         }

         for (m in 1 : n_chain)
         {
            param[[m]] = out[[m]]$param
            
            if (m == 1)
            {
               # save MCMC samples for the cold chain
               mcmc_samples$E0[ , , t]     = param[[1]]$E0
               mcmc_samples$D[ , , t]      = param[[1]]$D
               mcmc_samples$E1[ , , t]     = param[[1]]$E1
               mcmc_samples$U[ , , t]      = param[[1]]$U
               mcmc_samples$zeta[ , , t]   = param[[1]]$zeta
               mcmc_samples$eta[ , , t]    = param[[1]]$eta
               mcmc_samples$alpha0[ , , t] = param[[1]]$alpha0
               mcmc_samples$alpha1[ , , t] = param[[1]]$alpha1
               mcmc_samples$beta0[ , , t]  = param[[1]]$beta0
               mcmc_samples$beta1[ , , t]  = param[[1]]$beta1
               mcmc_samples$delta0[ , t]   = param[[1]]$delta0
               mcmc_samples$delta1[ , t]   = param[[1]]$delta1
               mcmc_samples$gamma0[ , t]   = param[[1]]$gamma0
               mcmc_samples$gamma1[ , t]   = param[[1]]$gamma1
               mcmc_samples$psi0[ , t]     = param[[1]]$psi0
               mcmc_samples$psi1[ , t]     = param[[1]]$psi1
               mcmc_samples$tau[ , t]      = param[[1]]$tau
               mcmc_samples$rho[ , t]      = param[[1]]$rho
               
               # save acceptance rates for the cold chain
               acceptance$zeta[ , , t]   = out[[1]]$acceptance$zeta
               acceptance$eta[ , , t]    = out[[1]]$acceptance$eta
               acceptance$alpha0[ , , t] = out[[1]]$acceptance$alpha0
               acceptance$alpha1[ , , t] = out[[1]]$acceptance$alpha1
               acceptance$beta0[ , , t]  = out[[1]]$acceptance$beta0
               acceptance$beta1[ , , t]  = out[[1]]$acceptance$beta1
               acceptance$delta0[ , t]   = out[[1]]$acceptance$delta0
               acceptance$delta1[ , t]   = out[[1]]$acceptance$delta1
               acceptance$gamma0[ , t]   = out[[1]]$acceptance$gamma0
               acceptance$gamma1[ , t]   = out[[1]]$acceptance$gamma1
               acceptance$psi0[ , t]     = out[[1]]$acceptance$psi0
               acceptance$psi1[ , t]     = out[[1]]$acceptance$psi1
            }
         }
      }
      else
      {
         # propose a swapping move
         # randomly choose chains to swap
         rand = sort(sample(1 : n_chain, 2))
         m1   = rand[1]
         m2   = rand[2]
         
         # calculate MH ratio for swapping
         ratio_swap = log_dDN(data0, data1, param[[m1]], priors, temperature[m2]) +
            log_dDN(data0, data1, param[[m2]], priors, temperature[m1]) -
            log_dDN(data0, data1, param[[m1]], priors, temperature[m1]) -
            log_dDN(data0, data1, param[[m2]], priors, temperature[m2])
         ratio_swap = exp(ratio_swap)
         
         # accept swapping
         if (is.nan(ratio_swap)) ratio_swap = 0
         if (runif(1) < min(1, ratio_swap))
         {
            #cat("swap the states of chains", m1, "and", m2, "\n")
            param_m1    = param[[m1]]
            param[[m1]] = param[[m2]]
            param[[m2]] = param_m1
         }
         
         # save MCMC samples for the cold chain
         mcmc_samples$E0[ , , t]     = param[[1]]$E0
         mcmc_samples$D[ , , t]      = param[[1]]$D
         mcmc_samples$E1[ , , t]     = param[[1]]$E1
         mcmc_samples$U[ , , t]      = param[[1]]$U
         mcmc_samples$zeta[ , , t]   = param[[1]]$zeta
         mcmc_samples$eta[ , , t]    = param[[1]]$eta
         mcmc_samples$alpha0[ , , t] = param[[1]]$alpha0
         mcmc_samples$alpha1[ , , t] = param[[1]]$alpha1
         mcmc_samples$beta0[ , , t]  = param[[1]]$beta0
         mcmc_samples$beta1[ , , t]  = param[[1]]$beta1
         mcmc_samples$delta0[ , t]   = param[[1]]$delta0
         mcmc_samples$delta1[ , t]   = param[[1]]$delta1
         mcmc_samples$gamma0[ , t]   = param[[1]]$gamma0
         mcmc_samples$gamma1[ , t]   = param[[1]]$gamma1
         mcmc_samples$psi0[ , t]     = param[[1]]$psi0
         mcmc_samples$psi1[ , t]     = param[[1]]$psi1
         mcmc_samples$tau[ , t]      = param[[1]]$tau
         mcmc_samples$rho[ , t]      = param[[1]]$rho
      }
      
      # print progress
      if (verbose)
      {
         if (t %% n_report == 0)
            cat("iter =", t, "\n")
      }
   }
   
   if (do_par) 
      stopCluster(cluster)
   
   # return a list of 
   # 1. MCMC samples (after burn-in period) for each parameter 
   # 2. Metropolis sampler acceptance rates for each parameter 
   results = list()
   results$samples = list(E0     = mcmc_samples$E0[ , , (n_burnin + 1) : n_sample],
                          D      = mcmc_samples$D[ , , (n_burnin + 1) : n_sample],
                          E1     = mcmc_samples$E1[ , , (n_burnin + 1) : n_sample],
                          U      = mcmc_samples$U[ , , (n_burnin + 1) : n_sample],
                          zeta   = mcmc_samples$zeta[ , , (n_burnin + 1) : n_sample],
                          eta    = mcmc_samples$eta[ , , (n_burnin + 1) : n_sample],
                          alpha0 = mcmc_samples$alpha0[ , , (n_burnin + 1) : n_sample],
                          alpha1 = mcmc_samples$alpha1[ , , (n_burnin + 1) : n_sample],
                          beta0  = mcmc_samples$beta0[ , , (n_burnin + 1) : n_sample],
                          beta1  = mcmc_samples$beta1[ , , (n_burnin + 1) : n_sample],
                          delta0 = mcmc_samples$delta0[ , (n_burnin + 1) : n_sample],
                          delta1 = mcmc_samples$delta1[ , (n_burnin + 1) : n_sample],
                          gamma0 = mcmc_samples$gamma0[ , (n_burnin + 1) : n_sample],
                          gamma1 = mcmc_samples$gamma1[ , (n_burnin + 1) : n_sample],
                          psi0   = mcmc_samples$psi0[ , (n_burnin + 1) : n_sample],
                          psi1   = mcmc_samples$psi1[ , (n_burnin + 1) : n_sample],
                          tau   = mcmc_samples$tau[ , (n_burnin + 1) : n_sample],
                          rho   = mcmc_samples$rho[ , (n_burnin + 1) : n_sample])
   results$acceptance = list(zeta   = 100 * mean(acceptance$zeta, na.rm = TRUE),
                             eta    = 100 * mean(acceptance$eta, na.rm = TRUE),
                             alpha0 = 100 * mean(acceptance$alpha0, na.rm = TRUE),
                             alpha1 = 100 * mean(acceptance$alpha1, na.rm = TRUE),
                             beta0  = 100 * mean(acceptance$beta0, na.rm = TRUE),
                             beta1  = 100 * mean(acceptance$beta1, na.rm = TRUE),
                             delta0 = 100 * mean(acceptance$delta0, na.rm = TRUE),
                             delta1 = 100 * mean(acceptance$delta1, na.rm = TRUE),
                             gamma0 = 100 * mean(acceptance$gamma0, na.rm = TRUE),
                             gamma1 = 100 * mean(acceptance$gamma1, na.rm = TRUE),
                             psi0   = 100 * mean(acceptance$psi0, na.rm = TRUE),
                             psi1   = 100 * mean(acceptance$psi1, na.rm = TRUE))
   class(results) = "2sampleDAG0"
   
   return(results)
}

#' Implementation of Gibbs step for two-sample DAG0
#'
#' @param param a list with each tag corresponding to a parameter name. 
#' Valid tags are 'E0', 'D', 'E1', 'U', 'zeta', 'eta' 'alpha0', 'alpha1', 'beta0', 'beta1', 'delta0', 'delta1', 'gamma0', 'gamma1', 'psi0', 'psi1', 'tau', and 'rho'. 
#' The value portion of each tag is the parameter values to be updated.
#' @param tuning a list with each tag corresponding to a parameter name. 
#' Valid tags are 'E', 'E_rev', 'U', 'zeta', 'eta', 'alpha0', 'alpha1', 'beta0', 'beta1', 'delta0', 'delta1', 'gamma0', 'gamma1', 'psi0', and 'psi1'. 
#' The value portion of each tag defines the standard deviations of Normal proposal distributions for Metropolis sampler for each parameter. 
#' The tag 'E' is for the birth or death update schemes for E = \{E_0, E_1\}, while 'E_rev' is for the reversal schemes.
#' The value portion of both tags should consist of the standard deviations of Gaussian proposals for alpha, beta, delta, gamma, zeta, and eta for an update of E.
#' Also, the value portion of the tag 'U' should include the standard deviations of Gaussian proposals for zeta and eta for an update of U.
#' @param priors a list with each tag corresponding to a parameter name. 
#' Valid tags are 'psi', 'tau', and 'rho'. The value portion of each tag defines the hyperparameters for two-sample DAG0.
#' @param x0 a matrix containing data for the control group.
#' @param x1 a matrix containing data for the case/treatment group.
#' @param temp a temperature.
#'
#' @return a list with the tags 'param' and 'acceptance'. 
#' The value portion of the tag 'param' gives parameter values after performing Gibbs step for two-sample DAG0.
#' The value portion of the tag 'acceptance' shows the acceptance indicators for each parameter.
#' @export 
#'
#' @examples
#' 
gibbs_2sampleDAG0 = function(param, tuning, priors, x0, x1, temp = 1)
{
   # the number of nodes
   p = ncol(x1)
   
   # parameters which will be updated
   alpha1   = param$alpha1
   alpha0   = param$alpha0
   beta1    = param$beta1
   beta0    = param$beta0
   delta1   = param$delta1
   delta0   = param$delta0
   gamma1   = param$gamma1
   gamma0   = param$gamma0
   psi1     = param$psi1
   psi0     = param$psi0
   eta      = param$zeta         # change label
   theta    = param$eta          # change label
   U        = param$U
   A1       = param$E1           # change label
   D        = param$D
   A0       = param$E0           # change label
   tau      = param$tau
   rho      = param$rho
   logitPi1 = param$logitPi1
   logitPi0 = param$logitPi0
   logMu1   = param$logMu1
   logMu0   = param$logMu0
   
   # initialize acceptance indicators
   acceptance = list()
   acceptance$alpha1 = acceptance$alpha0 = matrix(NA, p, p)
   acceptance$beta1  = acceptance$beta0  = matrix(NA, p, p)
   acceptance$delta1 = acceptance$delta0 = rep(NA, p)
   acceptance$gamma1 = acceptance$gamma0 = rep(NA, p)
   acceptance$psi1   = acceptance$psi0   = rep(NA, p)
   acceptance$zeta    = acceptance$eta  = matrix(NA, p, p)
   #acceptance$U = acceptance$A1 = acceptance$A0 = NA
   
   # update alpha1, alpha0
   for (j in 1 : p)
   {
      # logit(pi) of node j for each group
      logitPi1_old = logitPi1[ , j]
      logitPi0_old = logitPi0[ , j]
      
      for (k in 1 : p)
      {
         if (j == k) next
         
         # current alpha0 and alpha1
         alpha1_old = alpha1[j, k]
         alpha0_old = alpha0[j, k]
         
         if (A0[j, k] == 0)
         {
            if (A1[j, k] == 0)
            {
               next
            } else
            {
               # proposal
               alpha1_new = rnorm(1, mean = alpha1_old, sd = tuning$alpha1)
               logitPi1_new = logitPi1_old + x1[ , k] * (alpha1_new - alpha1_old)
               
               # calculate MH ratio
               llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1_old, logMu1[ , j], psi1[j])
               llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1[ , j], psi1[j])
               ratio_MH  = sum(llik1_new - llik1_old) - 0.5 * tau[1] * (alpha1_new * alpha1_new - alpha1_old * alpha1_old)
               ratio_MH  = exp(ratio_MH / temp)
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               acceptance$alpha1[j, k] = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  alpha1[j, k] = alpha1_new
                  logitPi1_old = logitPi1_new   # update logit(pi) of node j for group 1
                  acceptance$alpha1[j, k] = 1
               }
            }
         } else
         {
            if (A1[j, k] == 0)
            {
               # proposal
               alpha0_new = rnorm(1, mean = alpha0_old, sd = tuning$alpha0)
               logitPi0_new = logitPi0_old + x0[ , k] * (alpha0_new - alpha0_old)
               
               # calculate MH ratio
               llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0_old, logMu0[ , j], psi0[j])
               llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0[ , j], psi0[j])
               ratio_MH  = sum(llik0_new - llik0_old) - 0.5 * tau[1] * (alpha0_new * alpha0_new - alpha0_old * alpha0_old)
               ratio_MH  = exp(ratio_MH / temp)
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               acceptance$alpha0[j, k] = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  alpha0[j, k] = alpha0_new
                  logitPi0_old = logitPi0_new   # update logit(pi) of node j for group 0
                  acceptance$alpha0[j, k] = 1
               }
               
            } else
            {
               # proposal
               alpha0_new = rnorm(1, mean = alpha0_old, sd = tuning$alpha0)
               alpha1_new = alpha0_new + eta[j, k]
               logitPi0_new = logitPi0_old + x0[ , k] * (alpha0_new - alpha0_old)
               logitPi1_new = logitPi1_old + x1[ , k] * (alpha1_new - alpha1_old)
               
               # calculate MH ratio
               llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0_old, logMu0[ , j], psi0[j])
               llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0[ , j], psi0[j])
               llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1_old, logMu1[ , j], psi1[j])
               llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1[ , j], psi1[j])
               ratio_MH  = sum(llik0_new - llik0_old) + sum(llik1_new - llik1_old) - 
                  0.5 * tau[1] * (alpha0_new * alpha0_new - alpha0_old * alpha0_old)
               ratio_MH  = exp(ratio_MH / temp)
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               acceptance$alpha1[j, k] = acceptance$alpha0[j, k] = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  alpha0[j, k] = alpha0_new
                  alpha1[j, k] = alpha1_new
                  logitPi0_old = logitPi0_new   # update logit(pi) of node j for group 0
                  logitPi1_old = logitPi1_new   # update logit(pi) of node j for group 1
                  acceptance$alpha1[j, k] = acceptance$alpha0[j, k] = 1
               }
            }
         }
      }
      
      # save the updated logit(pi)'s of node j for both groups
      logitPi1[ , j] = logitPi1_old
      logitPi0[ , j] = logitPi0_old
   }
   
   # update beta1, beta0
   for (j in 1 : p)
   {
      # log(mu) of node j for each group
      logMu1_old = logMu1[ , j]
      logMu0_old = logMu0[ , j]
      
      for (k in 1 : p)
      {
         if (j == k) next
         
         # current alpha0 and alpha1
         beta1_old = beta1[j, k]
         beta0_old = beta0[j, k]
         
         if (A0[j, k] == 0)
         {
            if (A1[j, k] == 0)
            {
               next
            } else
            {
               # proposal
               beta1_new  = rnorm(1, mean = beta1_old, sd = tuning$beta1)
               logMu1_new = logMu1_old + x1[ , k] * (beta1_new - beta1_old)
               
               # calculate MH ratio
               llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1_old, psi1[j])
               llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1_new, psi1[j])
               ratio_MH  = sum(llik1_new - llik1_old) - 0.5 * tau[2] * (beta1_new * beta1_new - beta1_old * beta1_old)
               ratio_MH  = exp(ratio_MH / temp)
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               acceptance$beta1[j, k] = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  beta1[j, k] = beta1_new
                  logMu1_old  = logMu1_new   # update log(mu) of node j for group 1
                  acceptance$beta1[j, k] = 1
               }
            }
         } else
         {
            if (A1[j, k] == 0)
            {
               # proposal
               beta0_new  = rnorm(1, mean = beta0_old, sd = tuning$beta0)
               logMu0_new = logMu0_old + x0[ , k] * (beta0_new - beta0_old)
               
               # calculate MH ratio
               llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0_old, psi0[j])
               llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0_new, psi0[j])
               ratio_MH  = sum(llik0_new - llik0_old) - 0.5 * tau[2] * (beta0_new * beta0_new - beta0_old * beta0_old)
               ratio_MH  = exp(ratio_MH / temp)
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               acceptance$beta0[j, k] = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  beta0[j, k] = beta0_new
                  logMu0_old  = logMu0_new   # update log(mu) of node j for group 0
                  acceptance$beta0[j, k] = 1
               }
               
            } else
            {
               # proposal
               beta0_new  = rnorm(1, mean = beta0_old, sd = tuning$beta0)
               beta1_new  = beta0_new + theta[j, k]
               logMu0_new = logMu0_old + x0[ , k] * (beta0_new - beta0_old)
               logMu1_new = logMu1_old + x1[ , k] * (beta1_new - beta1_old)
               
               # calculate MH ratio
               llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0_old, psi0[j])
               llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0_new, psi0[j])
               llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1_old, psi1[j])
               llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1_new, psi1[j])
               ratio_MH  = sum(llik0_new - llik0_old) + sum(llik1_new - llik1_old) - 
                  0.5 * tau[2] * (beta0_new * beta0_new - beta0_old * beta0_old)
               ratio_MH  = exp(ratio_MH / temp)
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               acceptance$beta1[j, k] = acceptance$beta0[j, k] = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  beta0[j, k] = beta0_new
                  beta1[j, k] = beta1_new
                  logMu0_old  = logMu0_new   # update log(mu) of node j for group 0
                  logMu1_old  = logMu1_new   # update log(mu) of node j for group 1
                  acceptance$beta1[j, k] = acceptance$beta0[j, k] = 1
               }
            }
         }
      }
      
      # save the updated log(mu)'s of node j for both groups
      logMu1[ , j] = logMu1_old
      logMu0[ , j] = logMu0_old
   }
   
   # update delta1, delta0
   for (j in 1 : p)
   {
      # current delta1
      delta1_old = delta1[j]
      
      # proposal 
      delta1_new   = rnorm(1, mean = delta1_old, sd = tuning$delta1)
      logitPi1_new = logitPi1[ , j] + (delta1_new - delta1_old)
      
      # calculate MH ratio
      llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
      llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1[ , j], psi1[j])
      ratio_MH  = sum(llik1_new - llik1_old) - 0.5 * tau[3] * (delta1_new * delta1_new - delta1_old * delta1_old)
      ratio_MH  = exp(ratio_MH / temp)
      
      # accept the proposed values with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
      acceptance$delta1[j] = 0
      if (runif(1) < min(1, ratio_MH))
      {
         delta1[j]      = delta1_new
         logitPi1[ , j] = logitPi1_new     # update logit(pi) of node j for group 1
         acceptance$delta1[j] = 1
      }
   }
   
   for (j in 1 : p)
   {
      # current delta0
      delta0_old = delta0[j]
      
      # proposal 
      delta0_new   = rnorm(1, mean = delta0_old, sd = tuning$delta0)
      logitPi0_new = logitPi0[ , j] + (delta0_new - delta0_old)
      
      # calculate MH ratio
      llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
      llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0[ , j], psi0[j])
      ratio_MH  = sum(llik0_new - llik0_old) - 0.5 * tau[3] * (delta0_new * delta0_new - delta0_old * delta0_old)
      ratio_MH  = exp(ratio_MH / temp)
      
      # accept the proposed values with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
      acceptance$delta0[j] = 0
      if (runif(1) < min(1, ratio_MH))
      {
         delta0[j]      = delta0_new
         logitPi0[ , j] = logitPi0_new     # update logit(pi) of node j for group 0
         acceptance$delta0[j] = 1
      }
   }
   
   # update gamma1, gamma0
   for (j in 1 : p)
   {
      # current gamma1
      gamma1_old = gamma1[j]
      
      # proposal 
      gamma1_new = rnorm(1, mean = gamma1_old, sd = tuning$gamma1)
      logMu1_new = logMu1[ , j] + (gamma1_new - gamma1_old)
      
      # calculate MH ratio
      llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
      llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1_new, psi1[j])
      ratio_MH  = sum(llik1_new - llik1_old) - 0.5 * tau[4] * (gamma1_new * gamma1_new - gamma1_old * gamma1_old)
      ratio_MH  = exp(ratio_MH / temp)
      
      # accept the proposed values with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
      acceptance$gamma1[j] = 0
      if (runif(1) < min(1, ratio_MH))
      {
         gamma1[j]    = gamma1_new
         logMu1[ , j] = logMu1_new     # update log(mu) of node j for group 1
         acceptance$gamma1[j] = 1
      }
   }
   
   for (j in 1 : p)
   {
      # current gamma0
      gamma0_old = gamma0[j]
      
      # proposal 
      gamma0_new = rnorm(1, mean = gamma0_old, sd = tuning$gamma0)
      logMu0_new = logMu0[ , j] + (gamma0_new - gamma0_old)
      
      # calculate MH ratio
      llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
      llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0_new, psi0[j])
      ratio_MH  = sum(llik0_new - llik0_old) - 0.5 * tau[4] * (gamma0_new * gamma0_new - gamma0_old * gamma0_old)
      ratio_MH  = exp(ratio_MH / temp)
      
      # accept the proposed values with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
      acceptance$gamma0[j] = 0
      if (runif(1) < min(1, ratio_MH))
      {
         gamma0[j]    = gamma0_new
         logMu0[ , j] = logMu0_new     # update log(mu) of node j for group 0
         acceptance$gamma0[j] = 1
      }
   }

   # update psi1, psi0
   for (j in 1 : p)
   {
      # current log(psi1)
      logPsi1_old = log(psi1[j])
      
      # proposal
      logPsi1_new = rnorm(1, mean = logPsi1_old, sd = tuning$psi1)
      
      # calculate MH ratio
      llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], exp(logPsi1_old))
      llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], exp(logPsi1_new))
      ratio_MH  = sum(llik1_new - llik1_old) + priors$psi[1] * (logPsi1_new - logPsi1_old) - 
         priors$psi[2] * (exp(logPsi1_new) - exp(logPsi1_old))
      ratio_MH  = exp(ratio_MH / temp)
      
      # accept the proposed values with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
      acceptance$psi1[j] = 0
      if (runif(1) < min(1, ratio_MH))
      {
         psi1[j] = exp(logPsi1_new)   # transformed back to psi when updating
         acceptance$psi1[j] = 1
      }
   }
   
   for (j in 1 : p)
   {
      # current log(psi0)
      logPsi0_old = log(psi0[j])
      
      # proposal
      logPsi0_new = rnorm(1, mean = logPsi0_old, sd = tuning$psi0)
      
      # calculate MH ratio
      llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], exp(logPsi0_old))
      llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], exp(logPsi0_new))
      ratio_MH  = sum(llik0_new - llik0_old) + priors$psi[1] * (logPsi0_new - logPsi0_old) - 
         priors$psi[2] * (exp(logPsi0_new) - exp(logPsi0_old))
      ratio_MH  = exp(ratio_MH / temp)
      
      # accept the proposed values with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
      acceptance$psi0[j] = 0
      if (runif(1) < min(1, ratio_MH))
      {
         psi0[j] = exp(logPsi0_new)   # transformed back to psi when updating
         acceptance$psi0[j] = 1
      }
   }
   
   # update eta, theta
   id_upd = which(U == 1, arr.ind = TRUE)
   n_upd  = nrow(id_upd)
   if (n_upd > 0)
   {
      for (l in 1 : n_upd)
      {
         # index of updating theta
         j = id_upd[l, 1]
         k = id_upd[l, 2]
         
         # current alpha1 and eta
         alpha1_old = alpha1[j, k]
         eta_old    = eta[j, k]
         
         # proposal
         eta_new      = rnorm(1, mean = eta_old, sd = tuning$zeta)
         alpha1_new   = alpha0[j, k] + eta_new
         logitPi1_new = logitPi1[ , j] + x1[ , k] * (alpha1_new - alpha1_old)
         
         # calculate MH ratio
         llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
         llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1[ , j], psi1[j])
         ratio_MH  = sum(llik1_new - llik1_old) - 0.5 * tau[5] * (eta_new * eta_new - eta_old * eta_old)
         ratio_MH  = exp(ratio_MH / temp)
         
         # accept the proposed values with probability of min(1, MH ratio)
         if (is.nan(ratio_MH)) ratio_MH = 0
         acceptance$zeta[j, k] = 0
         if (runif(1) < min(1, ratio_MH))
         {
            eta[j, k]      = eta_new
            alpha1[j, k]   = alpha1_new
            logitPi1[ , j] = logitPi1_new   # update logit(pi) of node j for group 1
            acceptance$zeta[j, k] = 1
         } 
      }
      
      for (l in 1 : n_upd)
      {
         # index of updating theta
         j = id_upd[l, 1]
         k = id_upd[l, 2]
         
         # current beta1 and theta
         beta1_old  = beta1[j, k]
         theta_old  = theta[j, k]
         
         # proposal
         theta_new  = rnorm(1, mean = theta_old, sd = tuning$eta)
         beta1_new  = beta0[j, k] + theta_new
         logMu1_new   = logMu1[ , j] + x1[ , k] * (beta1_new - beta1_old)
         
         # calculate MH ratio
         llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
         llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1_new, psi1[j])
         ratio_MH  = sum(llik1_new - llik1_old) - 0.5 * tau[6] * (theta_new * theta_new - theta_old * theta_old)
         ratio_MH  = exp(ratio_MH / temp)
         
         # accept the proposed values with probability of min(1, MH ratio)
         if (is.nan(ratio_MH)) ratio_MH = 0
         acceptance$eta[j, k] = 0
         if (runif(1) < min(1, ratio_MH))
         {
            theta[j, k]  = theta_new
            beta1[j, k]  = beta1_new
            logMu1[ , j] = logMu1_new     # update log(mu) of node j for group 1
            acceptance$eta[j, k] = 1
         }
      }
   }
   
   # update U
   id_upd = which(A0 == 1 & A1 == 1, arr.ind = TRUE)
   n_upd  = nrow(id_upd)
   if (n_upd > 0)
   {
      #acceptance$U = 0
      
      for (l in 1 : n_upd)
      {
         # index of updating U
         j = id_upd[l, 1]
         k = id_upd[l, 2]
         
         # current alpha1, beta1, alpha0, beta0, eta, and theta
         alpha1_old = alpha1[j, k]
         beta1_old  = beta1[j, k]
         alpha0_old = alpha0[j, k]
         beta0_old  = beta0[j, k]
         eta_old    = eta[j, k]
         theta_old  = theta[j, k]
         
         if (U[j, k] == 0)
         {
            # proposal
            U_new      = 1
            eta_new    = rnorm(1, mean = 0, sd = tuning$U[1])
            theta_new  = rnorm(1, mean = 0, sd = tuning$U[2])
            alpha0_new = alpha0_old - 0.5 * eta_new
            beta0_new  = beta0_old - 0.5 * theta_new
            alpha1_new = alpha1_old + 0.5 * eta_new
            beta1_new  = beta1_old + 0.5 * theta_new
            logitPi0_new = logitPi0[ , j] + x0[ , k] * (alpha0_new - alpha0_old)
            logMu0_new   = logMu0[ , j] + x0[ , k] * (beta0_new - beta0_old)
            logitPi1_new = logitPi1[ , j] + x1[ , k] * (alpha1_new - alpha1_old)
            logMu1_new   = logMu1[ , j] + x1[ , k] * (beta1_new - beta1_old)
            
            # calculate MH ratio
            llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
            llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1_new, psi1[j])
            llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
            llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0_new, psi0[j])
            ratio_likel = sum(llik1_new - llik1_old) + sum(llik0_new - llik0_old)
            ratio_prior = -0.5 * tau[1] * (alpha0_new * alpha0_new - alpha0_old * alpha0_old) - 
               0.5 * tau[2] * (beta0_new * beta0_new - beta0_old * beta0_old) + 
               0.5 * log(tau[5]) - 0.5 * tau[5] * eta_new * eta_new + 
               0.5 * log(tau[6]) - 0.5 * tau[6] * theta_new * theta_new + 
               log(rho[1]) - log(1 - rho[1])
            ratio_prop  = log(tuning$U[1]) + 0.5 * (eta_new / tuning$U[1])^2 +
               log(tuning$U[2]) + 0.5 * (theta_new / tuning$U[2])^2
            ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop) 
         } else
         {
            # proposal
            U_new = eta_new = theta_new = 0
            alpha0_new = alpha1_new = 0.5 * (alpha0_old + alpha1_old)
            beta0_new  = beta1_new  = 0.5 * (beta0_old + beta1_old)
            logitPi0_new = logitPi0[ , j] + x0[ , k] * (alpha0_new - alpha0_old)
            logMu0_new   = logMu0[ , j] + x0[ , k] * (beta0_new - beta0_old)
            logitPi1_new = logitPi1[ , j] + x1[ , k] * (alpha1_new - alpha1_old)
            logMu1_new   = logMu1[ , j] + x1[ , k] * (beta1_new - beta1_old)
            
            # calculate MH ratio
            llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
            llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1_new, psi1[j])
            llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
            llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0_new, psi0[j])
            ratio_likel = sum(llik1_new - llik1_old) + sum(llik0_new - llik0_old)
            ratio_prior =  -0.5 * tau[1] * (alpha0_new * alpha0_new - alpha0_old * alpha0_old) - 
               0.5 * tau[2] * (beta0_new * beta0_new - beta0_old * beta0_old) + 
               log(1 - rho[1]) - 0.5 * log(tau[5]) + 0.5 * tau[5] * eta_old * eta_old - 
               0.5 * log(tau[6]) + 0.5 * tau[6] * theta_old * theta_old - log(rho[1])
            ratio_prop  = -log(tuning$U[1]) - 0.5 * (eta_old / tuning$U[1])^2 -
               log(tuning$U[2]) - 0.5 * (theta_old / tuning$U[2])^2
            ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop) 
         }
         
         # accept the proposed values with probability of min(1, MH ratio)
         if (is.nan(ratio_MH)) ratio_MH = 0
         if (runif(1) < min(1, ratio_MH))
         {
            U[j, k]      = U_new
            eta[j, k]    = eta_new
            theta[j, k]  = theta_new
            alpha0[j, k] = alpha0_new
            beta0[j, k]  = beta0_new
            alpha1[j, k] = alpha1_new
            beta1[j, k]  = beta1_new
            logitPi0[ , j] = logitPi0_new   # update logit(pi) of node j for group 0
            logMu0[ , j]   = logMu0_new     # update log(mu) of node j for group 0
            logitPi1[ , j] = logitPi1_new   # update logit(pi) of node j for group 1
            logMu1[ , j]   = logMu1_new     # update log(mu) of node j for group 1
            #acceptance$U = 1
         }
      }
   }
   
   if (runif(1) > 0.5)
   {
      # update edges for each group
      if (runif(1) > 0.5)
      {
         # update A1 & D, given A0 (propose birth or death)
         #acceptance$A1 = 0
         for (j in 1 : p)
         {
            for (k in 1 : p)
            {
               if (j == k) next
               
               # current alpha1, beta1, delta1, gamma1, eta, theta, and U
               alpha1_old = alpha1[j, k]
               beta1_old  = beta1[j, k]
               delta1_old = delta1[j]
               gamma1_old = gamma1[j]
               eta_old    = eta[j, k]
               theta_old  = theta[j, k]
               U_old      = U[j, k]
               
               if (A0[j, k] == 0)
               {
                  if (A1[j, k] == 0)
                  {
                     # proceed MH step unless addition of an edge makes a cycle
                     A1[j, k] = A1_new = 1
                     G1 = graph_from_adjacency_matrix(A1)
                     A1[j, k] = 0
                     if (!is.dag(G1)) next
                     
                     # proposal
                     D_new = 1
                     U_new = eta_new = theta_new = NA
                     alpha1_new = rnorm(1, mean = 0, sd = tuning$E[1])
                     beta1_new  = rnorm(1, mean = 0, sd = tuning$E[2])
                     delta1_new = rnorm(1, mean = delta1_old, sd = tuning$E[3])
                     gamma1_new = rnorm(1, mean = gamma1_old, sd = tuning$E[4])
                     logitPi1_new = logitPi1[ , j] + x1[ , k] * (alpha1_new - alpha1_old) + (delta1_new - delta1_old)
                     logMu1_new   = logMu1[ , j] + x1[ , k] * (beta1_new - beta1_old) + (gamma1_new - gamma1_old)
                     
                     # calculate MH ratio
                     llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
                     llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1_new, psi1[j])
                     ratio_likel = sum(llik1_new - llik1_old)
                     ratio_prior = 0.5 * log(tau[1]) - 0.5 * tau[1] * alpha1_new * alpha1_new + 
                        0.5 * log(tau[2]) - 0.5 * tau[2] * beta1_new * beta1_new + 
                        log(rho[2]) - log(1 - rho[2]) - 
                        0.5 * tau[3] * (delta1_new * delta1_new - delta1_old * delta1_old) - 
                        0.5 * tau[4] * (gamma1_new * gamma1_new - gamma1_old * gamma1_old)
                     ratio_prop  = log(tuning$E[1]) + 0.5 * (alpha1_new / tuning$E[1])^2 + 
                        log(tuning$E[2]) + 0.5 * (beta1_new / tuning$E[2])^2
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                     
                  } else
                  {
                     # proposal
                     A1_new = 0
                     D_new = 0 
                     U_new = eta_new = theta_new = NA
                     alpha1_new = beta1_new = 0 
                     delta1_new = rnorm(1, mean = delta1_old, sd = tuning$E[3])
                     gamma1_new = rnorm(1, mean = gamma1_old, sd = tuning$E[4])
                     logitPi1_new = logitPi1[ , j] + x1[ , k] * (alpha1_new - alpha1_old) + (delta1_new - delta1_old)
                     logMu1_new   = logMu1[ , j] + x1[ , k] * (beta1_new - beta1_old) + (gamma1_new - gamma1_old)
                     
                     # calculate MH ratio
                     llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
                     llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1_new, psi1[j])
                     ratio_likel = sum(llik1_new - llik1_old)
                     ratio_prior = log(1 - rho[2]) - 0.5 * log(tau[1]) + 0.5 * tau[1] * alpha1_old * alpha1_old - 
                        0.5 * log(tau[2]) + 0.5 * tau[2] * beta1_old * beta1_old - log(rho[2]) - 
                        0.5 * tau[3] * (delta1_new * delta1_new - delta1_old * delta1_old) - 
                        0.5 * tau[4] * (gamma1_new * gamma1_new - gamma1_old * gamma1_old)
                     ratio_prop = -log(tuning$E[1]) - 0.5 * (alpha1_old / tuning$E[1])^2 - 
                        log(tuning$E[2]) - 0.5 * (beta1_old / tuning$E[2])^2
                     ratio_MH   = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  }
               } else
               {
                  if (A1[j, k] == 0)
                  {
                     # proceed MH step unless addition of an edge makes a cycle
                     A1[j, k] = A1_new = 1
                     G1 = graph_from_adjacency_matrix(A1)
                     A1[j, k] = 0
                     if (!is.dag(G1)) next
                     
                     # proposal
                     D_new = 0
                     U_new = rbinom(1, size = 1, prob = 0.5)
                     eta_new    = U_new * rnorm(1, mean = 0, sd = tuning$E[5])
                     theta_new  = U_new * rnorm(1, mean = 0, sd = tuning$E[6])
                     alpha1_new = alpha0[j, k] + eta_new
                     beta1_new  = beta0[j, k] + theta_new
                     delta1_new = rnorm(1, mean = delta1_old, sd = tuning$E[3])
                     gamma1_new = rnorm(1, mean = gamma1_old, sd = tuning$E[4])
                     logitPi1_new = logitPi1[ , j] + x1[ , k] * (alpha1_new - alpha1_old) + (delta1_new - delta1_old)
                     logMu1_new   = logMu1[ , j] + x1[ , k] * (beta1_new - beta1_old) + (gamma1_new - gamma1_old)
                     
                     # calculate MH ratio
                     llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
                     llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1_new, psi1[j])
                     ratio_likel = sum(llik1_new - llik1_old) 
                     ratio_prior = log(1 - rho[1]) + log(1 - rho[2]) - log(rho[2]) + 
                        (U_new == 1) * (-log(1 - rho[1]) + 0.5 * log(tau[5]) - 0.5 * tau[5] * eta_new * eta_new +
                                           0.5 * log(tau[6]) - 0.5 * tau[6] * theta_new * theta_new + log(rho[1])) -
                        0.5 * tau[3] * (delta1_new * delta1_new - delta1_old * delta1_old) - 
                        0.5 * tau[4] * (gamma1_new * gamma1_new - gamma1_old * gamma1_old)
                     ratio_prop  = -log(0.5) + 
                        (U_new == 1) * (log(tuning$E[5]) + 0.5 * (eta_new / tuning$E[5])^2 +
                                           log(tuning$E[6]) + 0.5 * (theta_new / tuning$E[6])^2)
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  } else
                  {
                     # proposal
                     A1_new = 0
                     D_new  = 1
                     U_new = eta_new = theta_new = NA
                     alpha1_new = beta1_new = 0 
                     delta1_new = rnorm(1, mean = delta1_old, sd = tuning$E[3])
                     gamma1_new = rnorm(1, mean = gamma1_old, sd = tuning$E[4])
                     logitPi1_new = logitPi1[ , j] + x1[ , k] * (alpha1_new - alpha1_old) + (delta1_new - delta1_old)
                     logMu1_new   = logMu1[ , j] + x1[ , k] * (beta1_new - beta1_old) + (gamma1_new - gamma1_old)
                     
                     # calculate MH ratio
                     llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
                     llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1_new, psi1[j])
                     ratio_likel = sum(llik1_new - llik1_old) 
                     ratio_prior = log(rho[2]) - log(1 - rho[1]) - log(1 - rho[2]) + 
                        (U_old == 1) * (log(1 - rho[1]) - 0.5 * log(tau[5]) + 0.5 * tau[5] * eta_old * eta_old - 
                                           0.5 * log(tau[6]) + 0.5 * tau[6] * theta_old * theta_old - log(rho[1])) - 
                        0.5 * tau[3] * (delta1_new * delta1_new - delta1_old * delta1_old) - 
                        0.5 * tau[4] * (gamma1_new * gamma1_new - gamma1_old * gamma1_old)
                     ratio_prop  = log(0.5) + 
                        (U_old == 1) * (-log(tuning$E[5]) - 0.5 * (eta_old / tuning$E[5])^2 -
                                           log(tuning$E[6]) - 0.5 * (theta_old / tuning$E[6])^2 )
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  }
               }
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  A1[j, k] = A1_new
                  D[j, k]  = D_new
                  U[j, k]  = U_new
                  eta[j, k]    = eta_new
                  theta[j, k]  = theta_new
                  alpha1[j, k] = alpha1_new
                  beta1[j, k]  = beta1_new
                  delta1[j]    = delta1_new
                  gamma1[j]    = gamma1_new
                  logitPi1[ , j] = logitPi1_new   # update logit(pi) of node j for group 1
                  logMu1[ , j]   = logMu1_new     # update log(mu) of node j for group 1
                  #acceptance$A1 = 1
               }
            }
         }
         
         # update A0 & D, given A1 (propose birth or death)
         #acceptance$A0 = 0
         for (j in 1 : p)
         {
            for (k in 1 : p)
            {
               if (j == k) next
               
               # current alpha0, beta0, delta0, gamma0, eta, theta, and U
               alpha0_old = alpha0[j, k]
               beta0_old  = beta0[j, k]
               delta0_old = delta0[j]
               gamma0_old = gamma0[j]
               eta_old    = eta[j, k]
               theta_old  = theta[j, k]
               U_old      = U[j, k]
               
               if (A1[j, k] == 0)
               {
                  if (A0[j, k] == 0)
                  {
                     # proceed MH step unless our proposal makes a cycle
                     A0[j, k] = A0_new = 1
                     G0 = graph_from_adjacency_matrix(A0)
                     A0[j, k] = 0
                     if (!is.dag(G0)) next
                     
                     # proposal
                     D_new = 1
                     U_new = eta_new = theta_new = NA
                     alpha0_new = rnorm(1, mean = 0, sd = tuning$E[1])
                     beta0_new  = rnorm(1, mean = 0, sd = tuning$E[2])
                     delta0_new = rnorm(1, mean = delta0_old, sd = tuning$E[3])
                     gamma0_new = rnorm(1, mean = gamma0_old, sd = tuning$E[4])
                     logitPi0_new = logitPi0[ , j] + x0[ , k] * (alpha0_new - alpha0_old) + (delta0_new - delta0_old)
                     logMu0_new   = logMu0[ , j] + x0[ , k] * (beta0_new - beta0_old) + (gamma0_new - gamma0_old)
                     
                     # calculate MH ratio
                     llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
                     llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0_new, psi0[j])
                     ratio_likel = sum(llik0_new - llik0_old) 
                     ratio_prior = 0.5 * log(tau[1]) - 0.5 * tau[1] * alpha0_new * alpha0_new + 
                        0.5 * log(tau[2]) - 0.5 * tau[2] * beta0_new * beta0_new +
                        log(rho[2]) + log(rho[3]) - log(1 - rho[2]) - log(1 - rho[3]) - 
                        0.5 * tau[3] * (delta0_new * delta0_new - delta0_old * delta0_old) - 
                        0.5 * tau[4] * (gamma0_new * gamma0_new - gamma0_old * gamma0_old) 
                     ratio_prop  = log(tuning$E[1]) + 0.5 * (alpha0_new / tuning$E[1])^2 + 
                        log(tuning$E[2]) + 0.5 * (beta0_new / tuning$E[2])^2 
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  } else
                  {
                     # proposal
                     A0_new = 0
                     D_new  = 0
                     U_new  = eta_new = theta_new = NA
                     alpha0_new = beta0_new = 0
                     delta0_new = rnorm(1, mean = delta0_old, sd = tuning$E[3])
                     gamma0_new = rnorm(1, mean = gamma0_old, sd = tuning$E[4])
                     logitPi0_new = logitPi0[ , j] + x0[ , k] * (alpha0_new - alpha0_old) + (delta0_new - delta0_old)
                     logMu0_new   = logMu0[ , j] + x0[ , k] * (beta0_new - beta0_old) + (gamma0_new - gamma0_old)
                     
                     # calculate MH ratio
                     llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
                     llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0_new, psi0[j])
                     ratio_likel = sum(llik0_new - llik0_old) 
                     ratio_prior = log(1 - rho[2]) + log(1 - rho[3]) -
                        0.5 * log(tau[1]) + 0.5 * tau[1] * alpha0_old * alpha0_old -
                        0.5 * log(tau[2]) + 0.5 * tau[2] * beta0_old * beta0_old -
                        log(rho[2]) - log(rho[3]) - 
                        0.5 * tau[3] * (delta0_new * delta0_new - delta0_old * delta0_old) - 
                        0.5 * tau[4] * (gamma0_new * gamma0_new - gamma0_old * gamma0_old)
                     ratio_prop  = -log(tuning$E[1]) - 0.5 * (alpha0_old / tuning$E[1])^2 - 
                        log(tuning$E[2]) - 0.5 * (beta0_old / tuning$E[2])^2 
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  }
               } else
               {
                  if (A0[j, k] == 0)
                  {
                     # proceed MH step unless our proposal makes a cycle
                     A0[j, k] = A0_new = 1
                     G0 = graph_from_adjacency_matrix(A0)
                     A0[j, k] = 0
                     if (!is.dag(G0)) next
                     
                     # proposal
                     D_new = 0
                     U_new = rbinom(1, size = 1, prob = 0.5)
                     eta_new    = U_new * rnorm(1, mean = 0, sd = tuning$E[5])
                     theta_new  = U_new * rnorm(1, mean = 0, sd = tuning$E[6])
                     alpha0_new = alpha1[j, k] - eta_new
                     beta0_new  = beta1[j, k] - theta_new
                     delta0_new = rnorm(1, mean = delta0_old, sd = tuning$E[3])
                     gamma0_new = rnorm(1, mean = gamma0_old, sd = tuning$E[4])
                     logitPi0_new = logitPi0[ , j] + x0[ , k] * (alpha0_new - alpha0_old) + (delta0_new - delta0_old)
                     logMu0_new   = logMu0[ , j] + x0[ , k] * (beta0_new - beta0_old) + (gamma0_new - gamma0_old)
                     
                     # calculate MH ratio
                     llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
                     llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0_new, psi0[j])
                     ratio_likel = sum(llik0_new - llik0_old) 
                     ratio_prior -0.5 * tau[1] * (alpha0_new * alpha0_new - alpha1[j, k] * alpha1[j, k]) - 
                        0.5 * tau[2] * (beta0_new * beta0_new - beta1[j, k] * beta1[j, k]) + 
                        log(1 - rho[1]) + log(1 - rho[2]) + log(rho[3]) - log(rho[2]) - log(1 - rho[3]) + 
                        (U_new == 1) * (-log(1 - rho[1]) + 0.5 * log(tau[5]) - 0.5 * tau[5] * eta_new * eta_new + 
                                           0.5 * log(tau[6]) - 0.5 * tau[6] * theta_new * theta_new + log(rho[1])) - 
                        0.5 * tau[3] * (delta0_new * delta0_new - delta0_old * delta0_old) - 
                        0.5 * tau[4] * (gamma0_new * gamma0_new - gamma0_old * gamma0_old)
                     ratio_prop  = -log(0.5) + 
                        (U_new == 1) * (log(tuning$E[5]) + 0.5 * (eta_new / tuning$E[5])^2 +
                                           log(tuning$E[6]) + 0.5 * (theta_new / tuning$E[6])^2)
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  } else
                  {
                     # proposal
                     A0_new = 0
                     D_new  = 1
                     U_new  = eta_new = theta_new = NA
                     alpha0_new = beta0_new = 0
                     delta0_new = rnorm(1, mean = delta0_old, sd = tuning$E[3])
                     gamma0_new = rnorm(1, mean = gamma0_old, sd = tuning$E[4])
                     logitPi0_new = logitPi0[ , j] + x0[ , k] * (alpha0_new - alpha0_old) + (delta0_new - delta0_old)
                     logMu0_new   = logMu0[ , j] + x0[ , k] * (beta0_new - beta0_old) + (gamma0_new - gamma0_old)
                     
                     # calculate MH ratio
                     llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
                     llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0_new, psi0[j])
                     ratio_likel = sum(llik0_new - llik0_old) 
                     ratio_prior = -0.5 * tau[1] * (alpha1[j, k] * alpha1[j, k] - alpha0_old * alpha0_old) - 
                        0.5 * tau[2] * (beta1[j ,k] * beta1[j, k] - beta0_old * beta0_old) + 
                        log(rho[2]) + log(1 - rho[3]) - log(1 - rho[1]) - log(1 - rho[2]) - log(rho[3]) + 
                        (U_old == 1) * (log(1 - rho[1]) - 0.5 * log(tau[5]) + 0.5 * tau[5] * eta_old * eta_old - 
                                           0.5 * log(tau[6]) + 0.5 * tau[6] * theta_old * theta_old - log(rho[1])) - 
                        0.5 * tau[3] * (delta0_new * delta0_new - delta0_old * delta0_old) - 
                        0.5 * tau[4] * (gamma0_new * gamma0_new - gamma0_old * gamma0_old)
                     ratio_prop  = log(0.5) +
                        (U_old == 1) * (-log(tuning$E[5]) - 0.5 * (eta_old / tuning$E[5])^2 -
                                           log(tuning$E[6]) - 0.5 * (theta_old / tuning$E[6])^2)
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  }
               }
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  A0[j, k] = A0_new
                  D[j, k]  = D_new
                  U[j, k]  = U_new
                  eta[j, k]    = eta_new
                  theta[j, k]  = theta_new
                  alpha0[j, k] = alpha0_new
                  beta0[j, k]  = beta0_new
                  delta0[j]    = delta0_new
                  gamma0[j]    = gamma0_new
                  logitPi0[ , j] = logitPi0_new   # update logit(pi) of node j for group 0
                  logMu0[ , j]   = logMu0_new     # update log(mu) of node j for group 0
                  #acceptance$A0 = 1
               }
            }
         }
      } else
      {
         # update A1 & D, given A0 (propose reversal)
         id_rev = which(A1 == 1, arr.ind = TRUE)
         n_rev  = nrow(id_rev)
         if (n_rev > 0)
         {
            #acceptance$A1 = 0
            
            for (l in 1 : n_rev)
            {
               # index of an edge which will be reversed
               j = id_rev[l, 1]
               k = id_rev[l, 2]
               
               # propose reversal of the edge unless it makes a cycle
               A1[j, k] = A1_jk_new = 0
               A1[k, j] = A1_kj_new = 1
               G1 = graph_from_adjacency_matrix(A1)
               A1[j, k] = 1
               A1[k, j] = 0
               if (!is_dag(G1)) next
               
               # current alpha1, beta1, eta, theta and U corresponding to the reversal
               alpha1_jk_old = alpha1[j, k]
               alpha1_kj_old = alpha1[k, j]
               beta1_jk_old  = beta1[j, k]
               beta1_kj_old  = beta1[k, j]
               delta1_j_old  = delta1[j]
               delta1_k_old  = delta1[k]
               gamma1_j_old  = gamma1[j]
               gamma1_k_old  = gamma1[k]
               eta_jk_old    = eta[j, k]
               eta_kj_old    = eta[k, j]
               theta_jk_old  = theta[j, k]
               theta_kj_old  = theta[k, j]
               U_jk_old      = U[j, k]
               U_kj_old      = U[k, j]
               
               if (A0[j, k] == 0)
               {
                  if (A0[k, j] == 0)
                  {
                     # proposal
                     D_jk_new = 0
                     D_kj_new = 1
                     U_jk_new = eta_jk_new = theta_jk_new = NA
                     U_kj_new = eta_kj_new = theta_kj_new = NA
                     alpha1_jk_new = beta1_jk_new = 0
                     alpha1_kj_new = rnorm(1, mean = 0, sd = tuning$E_rev[1])
                     beta1_kj_new  = rnorm(1, mean = 0, sd = tuning$E_rev[2])
                     delta1_j_new  = rnorm(1, mean = delta1_j_old, sd = tuning$E_rev[3])
                     delta1_k_new  = rnorm(1, mean = delta1_k_old, sd = tuning$E_rev[3])
                     gamma1_j_new  = rnorm(1, mean = gamma1_j_old, sd = tuning$E_rev[4])
                     gamma1_k_new  = rnorm(1, mean = gamma1_k_old, sd = tuning$E_rev[4])
                     logitPi1_j_new = logitPi1[ , j] + x1[ , k] * (alpha1_jk_new - alpha1_jk_old) + (delta1_j_new - delta1_j_old)
                     logitPi1_k_new = logitPi1[ , k] + x1[ , j] * (alpha1_kj_new - alpha1_kj_old) + (delta1_k_new - delta1_k_old)
                     logMu1_j_new   = logMu1[ , j] + x1[ , k] * (beta1_jk_new - beta1_jk_old) + (gamma1_j_new - gamma1_j_old)
                     logMu1_k_new   = logMu1[ , k] + x1[ , j] * (beta1_kj_new - beta1_kj_old) + (gamma1_k_new - gamma1_k_old)
                     
                     # calculate MH ratio
                     llik1_j_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
                     llik1_k_old = llik_ZINBBN_j(x1[ , k], logitPi1[ , k], logMu1[ , k], psi1[k])
                     llik1_j_new = llik_ZINBBN_j(x1[ , j], logitPi1_j_new, logMu1_j_new, psi1[j])
                     llik1_k_new = llik_ZINBBN_j(x1[ , k], logitPi1_k_new, logMu1_k_new, psi1[k])
                     ratio_likel = sum(llik1_j_new - llik1_j_old) + sum(llik1_k_new - llik1_k_old)
                     ratio_prior = -0.5 * tau[1] * (alpha1_kj_new * alpha1_kj_new - alpha1_jk_old * alpha1_jk_old) - 
                        0.5 * tau[2] * (beta1_kj_new * beta1_kj_new - beta1_jk_old * beta1_jk_old) - 
                        0.5 * tau[3] * (delta1_j_new * delta1_j_new - delta1_j_old * delta1_j_old + 
                                           delta1_k_new * delta1_k_new - delta1_k_old * delta1_k_old) - 
                        0.5 * tau[4] * (gamma1_j_new * gamma1_j_new - gamma1_j_old * gamma1_j_old + 
                                           gamma1_k_new * gamma1_k_new - gamma1_k_old * gamma1_k_old)
                     ratio_prop  = -0.5 * ((alpha1_jk_old / tuning$E_rev[1])^2 - (alpha1_kj_new / tuning$E_rev[1])^2) - 
                        0.5 * ((beta1_jk_old / tuning$E_rev[2])^2 - (beta1_kj_new / tuning$E_rev[2])^2) 
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  } else
                  {
                     # proposal
                     D_jk_new = 0
                     D_kj_new = 0
                     U_jk_new = eta_jk_new = theta_jk_new = NA
                     U_kj_new = rbinom(1, size = 1, prob = 0.5)
                     eta_kj_new   = U_kj_new * rnorm(1, mean = 0, sd = tuning$E_rev[5])
                     theta_kj_new = U_kj_new * rnorm(1, mean = 0, sd = tuning$E_rev[6])
                     alpha1_jk_new = beta1_jk_new = 0
                     alpha1_kj_new = alpha0[k, j] + eta_kj_new
                     beta1_kj_new  = beta0[k, j] + theta_kj_new
                     delta1_j_new  = rnorm(1, mean = delta1_j_old, sd = tuning$E_rev[3])
                     delta1_k_new  = rnorm(1, mean = delta1_k_old, sd = tuning$E_rev[3])
                     gamma1_j_new  = rnorm(1, mean = gamma1_j_old, sd = tuning$E_rev[4])
                     gamma1_k_new  = rnorm(1, mean = gamma1_k_old, sd = tuning$E_rev[4])
                     logitPi1_j_new = logitPi1[ , j] + x1[ , k] * (alpha1_jk_new - alpha1_jk_old) + (delta1_j_new - delta1_j_old)
                     logitPi1_k_new = logitPi1[ , k] + x1[ , j] * (alpha1_kj_new - alpha1_kj_old) + (delta1_k_new - delta1_k_old)
                     logMu1_j_new   = logMu1[ , j] + x1[ , k] * (beta1_jk_new - beta1_jk_old) + (gamma1_j_new - gamma1_j_old)
                     logMu1_k_new   = logMu1[ , k] + x1[ , j] * (beta1_kj_new - beta1_kj_old) + (gamma1_k_new - gamma1_k_old)
                     
                     # calculate MH ratio
                     llik1_j_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
                     llik1_k_old = llik_ZINBBN_j(x1[ , k], logitPi1[ , k], logMu1[ , k], psi1[k])
                     llik1_j_new = llik_ZINBBN_j(x1[ , j], logitPi1_j_new, logMu1_j_new, psi1[j])
                     llik1_k_new = llik_ZINBBN_j(x1[ , k], logitPi1_k_new, logMu1_k_new, psi1[k])
                     ratio_likel = sum(llik1_j_new - llik1_j_old) + sum(llik1_k_new - llik1_k_old) 
                     ratio_prior = log(1 - rho[1]) + 2 * log(1 - rho[2]) -
                        0.5 * log(tau[1]) + 0.5 * tau[1] * alpha1_jk_old * alpha1_jk_old -
                        0.5 * log(tau[2]) + 0.5 * tau[2] * beta1_jk_old * beta1_jk_old - 2 * log(rho[2]) + 
                        (U_kj_new == 1) * (-log(1 - rho[1]) + 0.5 * log(tau[5]) - 0.5 * tau[5] * eta_kj_new * eta_kj_new +
                                              0.5 * log(tau[6]) - 0.5 * tau[6] * theta_kj_new * theta_kj_new + log(rho[1])) - 
                        0.5 * tau[3] * (delta1_j_new * delta1_j_new - delta1_j_old * delta1_j_old + 
                                           delta1_k_new * delta1_k_new - delta1_k_old * delta1_k_old) - 
                        0.5 * tau[4] * (gamma1_j_new * gamma1_j_new - gamma1_j_old * gamma1_j_old + 
                                           gamma1_k_new * gamma1_k_new - gamma1_k_old * gamma1_k_old)
                     ratio_prop  = -log(tuning$E_rev[1]) - 0.5 * (alpha1_jk_old / tuning$E_rev[1])^2 - 
                        log(tuning$E_rev[2]) - 0.5 * (beta1_jk_old / tuning$E_rev[2])^2 - log(0.5) + 
                        (U_kj_new == 1) * (log(tuning$E_rev[5]) + 0.5 * (eta_kj_new / tuning$E_rev[5])^2 +
                                              log(tuning$E_rev[6]) + 0.5 * (theta_kj_new / tuning$E_rev[6])^2)
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  }
               } else
               {
                  # proposal
                  D_jk_new = 1
                  D_kj_new = 1
                  U_jk_new = eta_jk_new = theta_jk_new = NA
                  U_kj_new = eta_kj_new = theta_kj_new = NA
                  alpha1_jk_new = beta1_jk_new = 0
                  alpha1_kj_new = rnorm(1, mean = 0, sd = tuning$E_rev[1])
                  beta1_kj_new  = rnorm(1, mean = 0, sd = tuning$E_rev[2])
                  delta1_j_new  = rnorm(1, mean = delta1_j_old, sd = tuning$E_rev[3])
                  delta1_k_new  = rnorm(1, mean = delta1_k_old, sd = tuning$E_rev[3])
                  gamma1_j_new  = rnorm(1, mean = gamma1_j_old, sd = tuning$E_rev[4])
                  gamma1_k_new  = rnorm(1, mean = gamma1_k_old, sd = tuning$E_rev[4])
                  logitPi1_j_new = logitPi1[ , j] + x1[ , k] * (alpha1_jk_new - alpha1_jk_old) + (delta1_j_new - delta1_j_old)
                  logitPi1_k_new = logitPi1[ , k] + x1[ , j] * (alpha1_kj_new - alpha1_kj_old) + (delta1_k_new - delta1_k_old)
                  logMu1_j_new   = logMu1[ , j] + x1[ , k] * (beta1_jk_new - beta1_jk_old) + (gamma1_j_new - gamma1_j_old)
                  logMu1_k_new   = logMu1[ , k] + x1[ , j] * (beta1_kj_new - beta1_kj_old) + (gamma1_k_new - gamma1_k_old)
                  
                  # calculate MH ratio
                  llik1_j_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
                  llik1_k_old = llik_ZINBBN_j(x1[ , k], logitPi1[ , k], logMu1[ , k], psi1[k])
                  llik1_j_new = llik_ZINBBN_j(x1[ , j], logitPi1_j_new, logMu1_j_new, psi1[j])
                  llik1_k_new = llik_ZINBBN_j(x1[ , k], logitPi1_k_new, logMu1_k_new, psi1[k])
                  ratio_likel = sum(llik1_j_new - llik1_j_old) + sum(llik1_k_new - llik1_k_old) 
                  ratio_prior = 0.5 * log(tau[1]) - 0.5 * tau[1] * alpha1_kj_new * alpha1_kj_new +
                     0.5 * log(tau[2]) - 0.5 * tau[2] * beta1_kj_new * beta1_kj_new + 2 * log(rho[2]) - 
                     log(1 - rho[1]) - 2 * log(1 - rho[2]) +
                     (U_jk_old == 1) * (log(1 - rho[1]) - 0.5 * log(tau[5]) + 0.5 * tau[5] * eta_jk_old * eta_jk_old -
                                           0.5 * log(tau[6]) + 0.5 * tau[6] * theta_jk_old * theta_jk_old - log(rho[1])) - 
                     0.5 * tau[3] * (delta1_j_new * delta1_j_new - delta1_j_old * delta1_j_old + 
                                        delta1_k_new * delta1_k_new - delta1_k_old * delta1_k_old) - 
                     0.5 * tau[4] * (gamma1_j_new * gamma1_j_new - gamma1_j_old * gamma1_j_old + 
                                        gamma1_k_new * gamma1_k_new - gamma1_k_old * gamma1_k_old)
                  ratio_prop  = log(0.5) + log(tuning$E_rev[1]) + 0.5 * (alpha1_kj_new / tuning$E_rev[1])^2 + 
                     log(tuning$E_rev[2]) + 0.5 * (beta1_kj_new / tuning$E_rev[2])^2 +
                     (U_jk_old == 1) * (-log(tuning$E_rev[5]) - 0.5 * (eta_jk_old / tuning$E_rev[5])^2 - 
                                           log(tuning$E_rev[6]) - 0.5 * (theta_jk_old / tuning$E_rev[6])^2) 
                  ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
               }
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  A1[j, k] = A1_jk_new
                  A1[k, j] = A1_kj_new
                  D[j, k]  = D_jk_new
                  D[k, j]  = D_kj_new
                  U[j, k]  = U_jk_new
                  U[k, j]  = U_kj_new
                  eta[j, k]    = eta_jk_new
                  eta[k, j]    = eta_kj_new
                  theta[j, k]  = theta_jk_new
                  theta[k, j]  = theta_kj_new
                  alpha1[j, k] = alpha1_jk_new
                  alpha1[k, j] = alpha1_kj_new
                  beta1[j, k]  = beta1_jk_new
                  beta1[k, j]  = beta1_kj_new
                  delta1[j]    = delta1_j_new
                  delta1[k]    = delta1_k_new
                  gamma1[j]    = gamma1_j_new
                  gamma1[k]    = gamma1_k_new
                  logitPi1[ , j] = logitPi1_j_new   # update logit(pi) of node j for group 1
                  logitPi1[ , k] = logitPi1_k_new   # update logit(pi) of node k for group 1
                  logMu1[ , j]   = logMu1_j_new     # update log(mu) of node j for group 1
                  logMu1[ , k]   = logMu1_k_new     # update log(mu) of node k for group 1
                  #acceptance$A1 = 1
               }
            }
         }
         
         # update A0 & D, given A1 (propose reversal)
         id_rev = which(A0 == 1, arr.ind = TRUE)
         n_rev  = nrow(id_rev)
         if (n_rev > 0)
         {
            #acceptance$A0 = 0
            for (l in 1 : n_rev)
            {
               # index of an edge which will be reversed
               j = id_rev[l, 1]
               k = id_rev[l, 2]
               
               # propose reversal of the edge unless it makes a cycle
               A0[j, k] = A0_jk_new = 0
               A0[k, j] = A0_kj_new = 1
               G0 = graph_from_adjacency_matrix(A0)
               A0[j, k] = 1
               A0[k, j] = 0
               if (!is_dag(G0)) next
               
               # current alpha0, beta0, eta, theta and U corresponding to the reversal
               alpha0_jk_old = alpha0[j, k]
               alpha0_kj_old = alpha0[k, j]
               beta0_jk_old  = beta0[j, k]
               beta0_kj_old  = beta0[k, j]
               delta0_j_old  = delta0[j]
               delta0_k_old  = delta0[k]
               gamma0_j_old  = gamma0[j]
               gamma0_k_old  = gamma0[k]
               eta_jk_old    = eta[j, k]
               eta_kj_old    = eta[k, j]
               theta_jk_old  = theta[j, k]
               theta_kj_old  = theta[k, j]
               U_jk_old      = U[j, k]
               U_kj_old      = U[k, j]
               
               if (A1[j, k] == 0)
               {
                  if (A1[k, j] == 0)
                  {
                     # proposal
                     D_jk_new = 0
                     D_kj_new = 1
                     U_jk_new = eta_jk_new = theta_jk_new = NA
                     U_kj_new = eta_kj_new = theta_kj_new = NA
                     alpha0_jk_new = beta0_jk_new = 0
                     alpha0_kj_new = rnorm(1, mean = 0, sd = tuning$E_rev[1])
                     beta0_kj_new  = rnorm(1, mean = 0, sd = tuning$E_rev[2])
                     delta0_j_new  = rnorm(1, mean = delta0_j_old, sd = tuning$E_rev[3])
                     delta0_k_new  = rnorm(1, mean = delta0_k_old, sd = tuning$E_rev[3])
                     gamma0_j_new  = rnorm(1, mean = gamma0_j_old, sd = tuning$E_rev[4])
                     gamma0_k_new  = rnorm(1, mean = gamma0_k_old, sd = tuning$E_rev[4])
                     logitPi0_j_new = logitPi0[ , j] + x0[ , k] * (alpha0_jk_new - alpha0_jk_old) + (delta0_j_new - delta0_j_old)
                     logitPi0_k_new = logitPi0[ , k] + x0[ , j] * (alpha0_kj_new - alpha0_kj_old) + (delta0_k_new - delta0_k_old)
                     logMu0_j_new   = logMu0[ , j] + x0[ , k] * (beta0_jk_new - beta0_jk_old) + (gamma0_j_new - gamma0_j_old)
                     logMu0_k_new   = logMu0[ , k] + x0[ , j] * (beta0_kj_new - beta0_kj_old) + (gamma0_k_new - gamma0_k_old)
                     
                     # calculate MH ratio
                     llik0_j_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
                     llik0_k_old = llik_ZINBBN_j(x0[ , k], logitPi0[ , k], logMu0[ , k], psi0[k])
                     llik0_j_new = llik_ZINBBN_j(x0[ , j], logitPi0_j_new, logMu0_j_new, psi0[j])
                     llik0_k_new = llik_ZINBBN_j(x0[ , k], logitPi0_k_new, logMu0_k_new, psi0[k])
                     ratio_likel = sum(llik0_j_new - llik0_j_old) + sum(llik0_k_new - llik0_k_old) 
                     ratio_prior = -0.5 * tau[1] * (alpha0_kj_new * alpha0_kj_new - alpha0_jk_old * alpha0_jk_old) -
                        0.5 * tau[2] * (beta0_kj_new * beta0_kj_new - beta0_jk_old * beta0_jk_old) - 
                        0.5 * tau[3] * (delta0_j_new * delta0_j_new - delta0_j_old * delta0_j_old + 
                                           delta0_k_new * delta0_k_new - delta0_k_old * delta0_k_old) - 
                        0.5 * tau[4] * (gamma0_j_new * gamma0_j_new - gamma0_j_old * gamma0_j_old + 
                                           gamma0_k_new * gamma0_k_new - gamma0_k_old * gamma0_k_old)
                     ratio_prop  = -0.5 * ((alpha0_jk_old / tuning$E_rev[1])^2 - (alpha0_kj_new / tuning$E_rev[1])^2) - 
                        0.5 * ((beta0_jk_old / tuning$E_rev[2])^2 - (beta0_kj_new / tuning$E_rev[2])^2)
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  } else
                  {
                     # proposal
                     D_jk_new = 0
                     D_kj_new = 0
                     U_jk_new = eta_jk_new = theta_jk_new = NA
                     U_kj_new = rbinom(1, size = 1, prob = 0.5)
                     eta_kj_new    = U_kj_new * rnorm(1, mean = 0 , sd = tuning$E_rev[5])
                     theta_kj_new  = U_kj_new * rnorm(1, mean = 0 , sd = tuning$E_rev[6])
                     alpha0_jk_new = beta0_jk_new = 0
                     alpha0_kj_new = alpha1[k, j] - eta_kj_new
                     beta0_kj_new  = beta1[k, j] - theta_kj_new
                     delta0_j_new  = rnorm(1, mean = delta0_j_old, sd = tuning$E_rev[3])
                     delta0_k_new  = rnorm(1, mean = delta0_k_old, sd = tuning$E_rev[3])
                     gamma0_j_new  = rnorm(1, mean = gamma0_j_old, sd = tuning$E_rev[4])
                     gamma0_k_new  = rnorm(1, mean = gamma0_k_old, sd = tuning$E_rev[4])
                     logitPi0_j_new = logitPi0[ , j] + x0[ , k] * (alpha0_jk_new - alpha0_jk_old) + (delta0_j_new - delta0_j_old)
                     logitPi0_k_new = logitPi0[ , k] + x0[ , j] * (alpha0_kj_new - alpha0_kj_old) + (delta0_k_new - delta0_k_old)
                     logMu0_j_new   = logMu0[ , j] + x0[ , k] * (beta0_jk_new - beta0_jk_old) + (gamma0_j_new - gamma0_j_old)
                     logMu0_k_new   = logMu0[ , k] + x0[ , j] * (beta0_kj_new - beta0_kj_old) + (gamma0_k_new - gamma0_k_old)
                     
                     # calculate MH ratio
                     llik0_j_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
                     llik0_k_old = llik_ZINBBN_j(x0[ , k], logitPi0[ , k], logMu0[ , k], psi0[k])
                     llik0_j_new = llik_ZINBBN_j(x0[ , j], logitPi0_j_new, logMu0_j_new, psi0[j])
                     llik0_k_new = llik_ZINBBN_j(x0[ , k], logitPi0_k_new, logMu0_k_new, psi0[k])
                     ratio_likel = sum(llik0_j_new - llik0_j_old) + sum(llik0_k_new - llik0_k_old) 
                     ratio_prior = -0.5 * tau[1] * (alpha0_kj_new * alpha0_kj_new - alpha0_jk_old * alpha0_jk_old) -
                        0.5 * tau[2] * (beta0_kj_new * beta0_kj_new - beta0_jk_old * beta0_jk_old) + 
                        log(1 - rho[1]) + 2 * log(1 - rho[2]) - 
                        0.5 * log(tau[1]) + 0.5 * tau[1] * alpha1[k, j] * alpha1[k, j] -
                        0.5 * log(tau[2]) + 0.5 * tau[2] * beta1[k, j] * beta1[k, j] - 2 * log(rho[2]) + 
                        (U_kj_new == 1) * (-log(1 - rho[1]) + 0.5 * log(tau[5]) - 0.5 * tau[5] * eta_kj_new * eta_kj_new +
                                              0.5 * log(tau[6]) - 0.5 * tau[6] * theta_kj_new * theta_kj_new + log(rho[1])) - 
                        0.5 * tau[3] * (delta0_j_new * delta0_j_new - delta0_j_old * delta0_j_old + 
                                           delta0_k_new * delta0_k_new - delta0_k_old * delta0_k_old) - 
                        0.5 * tau[4] * (gamma0_j_new * gamma0_j_new - gamma0_j_old * gamma0_j_old + 
                                           gamma0_k_new * gamma0_k_new - gamma0_k_old * gamma0_k_old)
                     ratio_prop  = -log(tuning$E_rev[1]) - 0.5 * (alpha0_jk_old / tuning$E_rev[1])^2 -
                        log(tuning$E_rev[2]) - 0.5 * (beta0_jk_old / tuning$E_rev[2])^2 - log(0.5) + 
                        (U_kj_new == 1) * (log(tuning$E_rev[5]) + 0.5 * (eta_kj_new / tuning$E_rev[5])^2 +
                                              log(tuning$E_rev[6]) + 0.5 * (theta_kj_new / tuning$E_rev[6])^2)
                     ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
                  }
               } else
               {
                  # proposal
                  D_jk_new = 1
                  D_kj_new = 1
                  U_jk_new = eta_jk_new = theta_jk_new = NA
                  U_kj_new = eta_kj_new = theta_kj_new = NA
                  alpha0_jk_new = beta0_jk_new = 0
                  alpha0_kj_new = rnorm(1, mean = 0, sd = tuning$E_rev[1])
                  beta0_kj_new  = rnorm(1, mean = 0, sd = tuning$E_rev[2])
                  delta0_j_new  = rnorm(1, mean = delta0_j_old, sd = tuning$E_rev[3])
                  delta0_k_new  = rnorm(1, mean = delta0_k_old, sd = tuning$E_rev[3])
                  gamma0_j_new  = rnorm(1, mean = gamma0_j_old, sd = tuning$E_rev[4])
                  gamma0_k_new  = rnorm(1, mean = gamma0_k_old, sd = tuning$E_rev[4])
                  logitPi0_j_new = logitPi0[ , j] + x0[ , k] * (alpha0_jk_new - alpha0_jk_old) + (delta0_j_new - delta0_j_old)
                  logitPi0_k_new = logitPi0[ , k] + x0[ , j] * (alpha0_kj_new - alpha0_kj_old) + (delta0_k_new - delta0_k_old)
                  logMu0_j_new   = logMu0[ , j] + x0[ , k] * (beta0_jk_new - beta0_jk_old) + (gamma0_j_new - gamma0_j_old)
                  logMu0_k_new   = logMu0[ , k] + x0[ , j] * (beta0_kj_new - beta0_kj_old) + (gamma0_k_new - gamma0_k_old)
                  
                  # calculate MH ratio
                  llik0_j_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
                  llik0_k_old = llik_ZINBBN_j(x0[ , k], logitPi0[ , k], logMu0[ , k], psi0[k])
                  llik0_j_new = llik_ZINBBN_j(x0[ , j], logitPi0_j_new, logMu0_j_new, psi0[j])
                  llik0_k_new = llik_ZINBBN_j(x0[ , k], logitPi0_k_new, logMu0_k_new, psi0[k])
                  ratio_likel = sum(llik0_j_new - llik0_j_old) + sum(llik0_k_new - llik0_k_old) 
                  ratio_prior = -0.5 * tau[1] * (alpha0_kj_new * alpha0_kj_new - alpha0_jk_old * alpha0_jk_old) -
                     0.5 * tau[2] * (beta0_kj_new * beta0_kj_new - beta0_jk_old * beta0_jk_old) +
                     0.5 * log(tau[1]) - 0.5 * tau[1] * alpha1[j, k] * alpha1[j, k] + 
                     0.5 * log(tau[2]) - 0.5 * tau[2] * beta1[j, k] * beta1[j, k] + 2 * log(rho[2]) -
                     log(1 - rho[1]) - 2 * log(1 - rho[2]) + 
                     (U_jk_old == 1) * (log(1 - rho[1]) - 0.5 * log(tau[5]) + 0.5 * tau[5] * eta_jk_old * eta_jk_old -
                                           0.5 * log(tau[6]) + 0.5 * tau[6] * theta_jk_old * theta_jk_old - log(rho[1])) - 
                     0.5 * tau[3] * (delta0_j_new * delta0_j_new - delta0_j_old * delta0_j_old + 
                                        delta0_k_new * delta0_k_new - delta0_k_old * delta0_k_old) - 
                     0.5 * tau[4] * (gamma0_j_new * gamma0_j_new - gamma0_j_old * gamma0_j_old + 
                                        gamma0_k_new * gamma0_k_new - gamma0_k_old * gamma0_k_old)
                  ratio_prop  = log(0.5) + log(tuning$E_rev[1]) + 0.5 * (alpha0_kj_new / tuning$E_rev[1])^2 +
                     log(tuning$E_rev[2]) + 0.5 * (beta0_kj_new / tuning$E_rev[2])^2 + 
                     (U_jk_old == 1) * (-log(tuning$E_rev[5]) - 0.5 * (eta_jk_old / tuning$E_rev[5])^2 - 
                                           log(tuning$E_rev[6]) - 0.5 * (theta_jk_old / tuning$E_rev[6])^2)
                  ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
               }
               
               # accept the proposed values with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  A0[j, k] = A0_jk_new
                  A0[k, j] = A0_kj_new
                  D[j, k]  = D_jk_new
                  D[k, j]  = D_kj_new
                  U[j, k]  = U_jk_new
                  U[k, j]  = U_kj_new
                  eta[j, k]    = eta_jk_new
                  eta[k, j]    = eta_kj_new
                  theta[j, k]  = theta_jk_new
                  theta[k, j]  = theta_kj_new
                  alpha0[j, k] = alpha0_jk_new
                  alpha0[k, j] = alpha0_kj_new
                  beta0[j, k]  = beta0_jk_new
                  beta0[k, j]  = beta0_kj_new
                  delta0[j]    = delta0_j_new
                  delta0[k]    = delta0_k_new
                  gamma0[j]    = gamma0_j_new
                  gamma0[k]    = gamma0_k_new
                  logitPi0[ , j] = logitPi0_j_new   # update logit(pi) of node j for group 0
                  logitPi0[ , k] = logitPi0_k_new   # update logit(pi) of node k for group 0
                  logMu0[ , j]   = logMu0_j_new     # update log(mu) of node j for group 0
                  logMu0[ , k]   = logMu0_k_new     # update log(mu) of node k for group 0
                  #acceptance$A0 = 1
               }
            }
         }
      }
   } else
   {
      # update shared edges for both groups
      if (runif(1) > 0.5)
      {
         # update A1 & A0, given D (propose birth or death of shared edges)
         for (j in 1 : p)
         {
            for (k in 1 : p)
            {
               if (j == k) next
               if (A1[j, k] != A0[j, k]) next
               
               # current 
               alpha1_old = alpha1[j, k]
               alpha0_old = alpha0[j, k]
               beta1_old  = beta1[j, k]
               beta0_old  = beta0[j, k]
               delta1_old = delta1[j]
               delta0_old = delta0[j]
               gamma1_old = gamma1[j]
               gamma0_old = gamma0[j]
               eta_old    = eta[j, k]
               theta_old  = theta[j, k] 
               U_old      = U[j, k]
               
               if (A0[j, k] == 0)
               {
                  # proceed MH step unless addition of a shared edge makes a cycle in A0 or A1 
                  A0[j, k] = A1[j, k] = A0_new = A1_new = 1
                  G0 = graph_from_adjacency_matrix(A0)
                  G1 = graph_from_adjacency_matrix(A1)
                  A0[j, k] = A1[j, k] = 0
                  if (!is.dag(G0) | !is.dag(G1)) next
                  
                  # proposal
                  U_new      = rbinom(1, size = 1, prob = 0.5)
                  eta_new    = U_new * rnorm(1, mean = 0, sd = tuning$E[5])
                  theta_new  = U_new * rnorm(1, mean = 0, sd = tuning$E[6])
                  alpha0_new = rnorm(1, mean = 0, sd = tuning$E[1])
                  alpha1_new = alpha0_new + eta_new
                  beta0_new  = rnorm(1, mean = 0, sd = tuning$E[2])
                  beta1_new  = beta0_new + theta_new
                  delta0_new = rnorm(1, mean = delta0_old, sd = tuning$E[3])
                  delta1_new = rnorm(1, mean = delta1_old, sd = tuning$E[3])
                  gamma0_new = rnorm(1, mean = gamma0_old, sd = tuning$E[4])
                  gamma1_new = rnorm(1, mean = gamma1_old, sd = tuning$E[4])
                  logitPi0_new = logitPi0[ , j] + x0[ , k] * (alpha0_new - alpha0_old) + (delta0_new - delta0_old)
                  logitPi1_new = logitPi1[ , j] + x1[ , k] * (alpha1_new - alpha1_old) + (delta1_new - delta1_old)
                  logMu0_new   = logMu0[ , j] + x0[ , k] * (beta0_new - beta0_old) + (gamma0_new - gamma0_old)
                  logMu1_new   = logMu1[ , j] + x1[ , k] * (beta1_new - beta1_old) + (gamma1_new - gamma1_old)
                  
                  # calculate MH ratio
                  llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
                  llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
                  llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0_new, psi0[j])
                  llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1_new, psi1[j])
                  ratio_likel = sum(llik1_new - llik1_old) + sum(llik0_new - llik0_old)
                  ratio_prior = 0.5 * log(tau[1]) - 0.5 * tau[1] * alpha0_new * alpha0_new + 
                     0.5 * log(tau[2]) - 0.5 * tau[2] * beta0_new * beta0_new + 
                     log(1 - rho[1]) + log(rho[3]) - log(1 - rho[3]) +
                     (U_new == 1) * (-log(1 - rho[1]) + 0.5 * log(tau[5]) - 0.5 * tau[5] * eta_new * eta_new + 
                                        0.5 * log(tau[6]) - 0.5 * tau[6] * theta_new * theta_new + log(rho[1])) - 
                     0.5 * tau[3] * (delta1_new * delta1_new - delta1_old * delta1_old + 
                                        delta0_new * delta0_new - delta0_old * delta0_old) - 
                     0.5 * tau[4] * (gamma1_new * gamma1_new - gamma1_old * gamma1_old +
                                        gamma0_new * gamma0_new - gamma0_old * gamma0_old)
                  ratio_prop  = -log(0.5) + log(tuning$E[1]) + 0.5 * (alpha0_new / tuning$E[1])^2 +
                     log(tuning$E[2]) + 0.5 * (beta0_new / tuning$E[2])^2 + 
                     (U_new == 1) * (log(tuning$E[5]) + 0.5 * (eta_new / tuning$E[5])^2 + 
                                        log(tuning$E[6]) + 0.5 * (theta_new / tuning$E[6])^2)
                  ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
               } else
               {
                  # proposal
                  A0_new = A1_new = 0
                  U_new = eta_new = theta_new = NA
                  alpha0_new = alpha1_new = 0
                  beta0_new  = beta1_new  = 0
                  delta0_new = rnorm(1, mean = delta0_old, sd = tuning$E[3])
                  delta1_new = rnorm(1, mean = delta1_old, sd = tuning$E[3])
                  gamma0_new = rnorm(1, mean = gamma0_old, sd = tuning$E[4])
                  gamma1_new = rnorm(1, mean = gamma1_old, sd = tuning$E[4])
                  logitPi0_new = logitPi0[ , j] + x0[ , k] * (alpha0_new - alpha0_old) + (delta0_new - delta0_old)
                  logitPi1_new = logitPi1[ , j] + x1[ , k] * (alpha1_new - alpha1_old) + (delta1_new - delta1_old)
                  logMu0_new   = logMu0[ , j] + x0[ , k] * (beta0_new - beta0_old) + (gamma0_new - gamma0_old)
                  logMu1_new   = logMu1[ , j] + x1[ , k] * (beta1_new - beta1_old) + (gamma1_new - gamma1_old)
                  
                  # calculate MH ratio
                  llik0_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
                  llik1_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
                  llik0_new = llik_ZINBBN_j(x0[ , j], logitPi0_new, logMu0_new, psi0[j])
                  llik1_new = llik_ZINBBN_j(x1[ , j], logitPi1_new, logMu1_new, psi1[j])
                  ratio_likel = sum(llik1_new - llik1_old) + sum(llik0_new - llik0_old)
                  ratio_prior = -0.5 * log(tau[1]) + 0.5 * tau[1] * alpha0_old * alpha0_old - 
                     0.5 * log(tau[2]) + 0.5 * tau[2] * beta0_old * beta0_old - 
                     log(1 - rho[1]) - log(rho[3]) + log(1 - rho[3]) +
                     (U_old == 1) * (log(1 - rho[1]) - 0.5 * log(tau[5]) + 0.5 * tau[5] * eta_old * eta_old - 
                                        0.5 * log(tau[6]) + 0.5 * tau[6] * theta_old * theta_old - log(rho[1])) - 
                     0.5 * tau[3] * (delta1_new * delta1_new - delta1_old * delta1_old + 
                                        delta0_new * delta0_new - delta0_old * delta0_old) - 
                     0.5 * tau[4] * (gamma1_new * gamma1_new - gamma1_old * gamma1_old +
                                        gamma0_new * gamma0_new - gamma0_old * gamma0_old)
                  ratio_prop  = log(0.5) - log(tuning$E[1]) - 0.5 * (alpha0_old / tuning$E[1])^2 -
                     log(tuning$E[2]) - 0.5 * (beta0_old / tuning$E[2])^2 +
                     (U_old == 1) * (-log(tuning$E[5]) - 0.5 * (eta_old / tuning$E[5])^2 -
                                        log(tuning$E[6]) - 0.5 * (theta_old / tuning$E[6])^2)
                  ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
               }
               
               # accept proposal with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  A0[j, k] = A0_new
                  A1[j, k] = A1_new
                  U[j, k]  = U_new
                  eta[j, k]    = eta_new
                  theta[j, k]  = theta_new
                  alpha0[j, k] = alpha0_new
                  alpha1[j, k] = alpha1_new
                  beta0[j, k]  = beta0_new
                  beta1[j, k]  = beta1_new
                  delta0[j]    = delta0_new
                  delta1[j]    = delta1_new
                  gamma0[j]    = gamma0_new
                  gamma1[j]    = gamma1_new
                  logitPi0[ , j] = logitPi0_new
                  logitPi1[ , j] = logitPi1_new
                  logMu0[ , j]   = logMu0_new
                  logMu1[ , j]   = logMu1_new
               }
            }
         }
      } else
      {
         # update A1 & A0, given D (propose reversal of shared edges)
         id_rev = which(A0 * A1 == 1, arr.ind = TRUE)
         n_rev  = nrow(id_rev)
         if (n_rev > 0)
         {
            for (l in 1 : n_rev)
            {
               # index of an edge which will be reversed
               j = id_rev[l, 1]
               k = id_rev[l, 2]
               
               # propose reversal of a shared edge unless it makes a cycle in A0 or A1
               A0[j, k] = A1[j, k] = A0_jk_new = A1_jk_new = 0
               A0[k, j] = A1[k, j] = A0_kj_new = A1_kj_new = 1
               G0 = graph_from_adjacency_matrix(A0)
               G1 = graph_from_adjacency_matrix(A1)
               A0[j, k] = A1[j, k] = 1
               A0[k, j] = A1[k, j] = 0
               if (!is.dag(G0) | !is.dag(G1)) next
               
               # current
               alpha1_jk_old = alpha1[j, k]
               alpha1_kj_old = alpha1[k, j]
               alpha0_jk_old = alpha0[j, k]
               alpha0_kj_old = alpha0[k, j]
               beta1_jk_old  = beta1[j, k]
               beta1_kj_old  = beta1[k, j]
               beta0_jk_old  = beta0[j, k]
               beta0_kj_old  = beta0[k, j]
               delta1_j_old  = delta1[j]
               delta1_k_old  = delta1[k]
               delta0_j_old  = delta0[j]
               delta0_k_old  = delta0[k]
               gamma1_j_old  = gamma1[j]
               gamma1_k_old  = gamma1[k]
               gamma0_j_old  = gamma0[j]
               gamma0_k_old  = gamma0[k]
               eta_jk_old    = eta[j, k]
               eta_kj_old    = eta[k, j]
               theta_jk_old  = theta[j, k]
               theta_kj_old  = theta[k, j]
               U_jk_old      = U[j, k]
               U_kj_old      = U[k, j]
               
               # proposal
               U_jk_new = eta_jk_new = theta_jk_new = NA
               U_kj_new = rbinom(1, size = 1, prob = 0.5)
               eta_kj_new    = U_kj_new * rnorm(1, mean = 0, sd = tuning$E_rev[5])
               theta_kj_new  = U_kj_new * rnorm(1, mean = 0, sd = tuning$E_rev[6])
               alpha0_jk_new = alpha1_jk_new = 0
               alpha0_kj_new = rnorm(1, mean = 0, sd = tuning$E_rev[1])
               alpha1_kj_new = alpha0_kj_new + eta_kj_new
               beta0_jk_new  = beta1_jk_new  = 0
               beta0_kj_new  = rnorm(1, mean = 0, sd = tuning$E_rev[2])
               beta1_kj_new  = beta0_kj_new + theta_kj_new
               delta0_j_new  = rnorm(1, mean = delta0_j_old, sd = tuning$E_rev[3])
               delta0_k_new  = rnorm(1, mean = delta0_k_old, sd = tuning$E_rev[3])
               delta1_j_new  = rnorm(1, mean = delta1_j_old, sd = tuning$E_rev[3])
               delta1_k_new  = rnorm(1, mean = delta1_k_old, sd = tuning$E_rev[3])
               gamma0_j_new  = rnorm(1, mean = gamma0_j_old, sd = tuning$E_rev[4])
               gamma0_k_new  = rnorm(1, mean = gamma0_k_old, sd = tuning$E_rev[4])
               gamma1_j_new  = rnorm(1, mean = gamma1_j_old, sd = tuning$E_rev[4])
               gamma1_k_new  = rnorm(1, mean = gamma1_k_old, sd = tuning$E_rev[4])
               logitPi0_j_new = logitPi0[ , j] + x0[ , k] * (alpha0_jk_new - alpha0_jk_old) + (delta0_j_new - delta0_j_old)
               logitPi0_k_new = logitPi0[ , k] + x0[ , j] * (alpha0_kj_new - alpha0_kj_old) + (delta0_k_new - delta0_k_old)
               logitPi1_j_new = logitPi1[ , j] + x1[ , k] * (alpha1_jk_new - alpha1_jk_old) + (delta1_j_new - delta1_j_old)
               logitPi1_k_new = logitPi1[ , k] + x1[ , j] * (alpha1_kj_new - alpha1_kj_old) + (delta1_k_new - delta1_k_old)
               logMu0_j_new   = logMu0[ , j] + x0[ , k] * (beta0_jk_new - beta0_jk_old) + (gamma0_j_new - gamma0_j_old)
               logMu0_k_new   = logMu0[ , k] + x0[ , j] * (beta0_kj_new - beta0_kj_old) + (gamma0_k_new - gamma0_k_old)
               logMu1_j_new   = logMu1[ , j] + x1[ , k] * (beta1_jk_new - beta1_jk_old) + (gamma1_j_new - gamma1_j_old)
               logMu1_k_new   = logMu1[ , k] + x1[ , j] * (beta1_kj_new - beta1_kj_old) + (gamma1_k_new - gamma1_k_old)

               # calculate MH ratio
               llik0_j_old = llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j])
               llik0_k_old = llik_ZINBBN_j(x0[ , k], logitPi0[ , k], logMu0[ , k], psi0[k])
               llik1_j_old = llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])
               llik1_k_old = llik_ZINBBN_j(x1[ , k], logitPi1[ , k], logMu1[ , k], psi1[k])
               llik0_j_new = llik_ZINBBN_j(x0[ , j], logitPi0_j_new, logMu0_j_new, psi0[j])
               llik0_k_new = llik_ZINBBN_j(x0[ , k], logitPi0_k_new, logMu0_k_new, psi0[k])
               llik1_j_new = llik_ZINBBN_j(x1[ , j], logitPi1_j_new, logMu1_j_new, psi1[j])
               llik1_k_new = llik_ZINBBN_j(x1[ , k], logitPi1_k_new, logMu1_k_new, psi1[k])
               ratio_likel = sum(llik1_j_new - llik1_j_old) + sum(llik1_k_new - llik1_k_old) +
                  sum(llik0_j_new - llik0_j_old) + sum(llik0_k_new - llik0_k_old)
               ratio_prior = -0.5 * tau[1] * (alpha0_kj_new *alpha0_kj_new - alpha0_jk_old * alpha0_jk_old) - 
                  0.5 * tau[2] * (beta0_kj_new * beta0_kj_new - beta0_jk_old * beta0_jk_old) +
                  (U_kj_new == 1) * (-log(1 - rho[1]) + 0.5 * log(tau[5]) - 0.5 * tau[5] * eta_kj_new * eta_kj_new + 
                                        0.5 * log(tau[6]) - 0.5 * tau[6] * theta_kj_new * theta_kj_new + log(rho[1])) -
                  (U_jk_old == 1) * (-log(1 - rho[1]) + 0.5 * log(tau[5]) - 0.5 * tau[5] * eta_jk_old * eta_jk_old + 
                                        0.5 * log(tau[6]) - 0.5 * tau[6] * theta_jk_old * theta_jk_old + log(rho[1])) - 
                  0.5 * tau[3] * (delta1_j_new * delta1_j_new - delta1_j_old * delta1_j_old +
                                     delta1_k_new * delta1_k_new - delta1_k_old * delta1_k_old +
                                     delta0_j_new * delta0_j_new - delta0_j_old * delta0_j_old + 
                                     delta0_k_new * delta0_k_new - delta0_k_old * delta0_k_old) - 
                  0.5 * tau[4] * (gamma1_j_new * gamma1_j_new - gamma1_j_old * gamma1_j_old +
                                     gamma1_k_new * gamma1_k_new - gamma1_k_old * gamma1_k_old +
                                     gamma0_j_new * gamma0_j_new - gamma0_j_old * gamma0_j_old + 
                                     gamma0_k_new * gamma0_k_new - gamma0_k_old * gamma0_k_old)
               ratio_prop  = -0.5 * (alpha0_jk_old / tuning$E_rev[1])^2 + 0.5 * (alpha0_kj_new / tuning$E_rev[1])^2 -
                  0.5 * (beta0_jk_old / tuning$E_rev[2])^2 + 0.5 * (beta0_kj_new / tuning$E_rev[2])^2 + 
                  (U_jk_old == 1) * (-log(tuning$E_rev[5]) - 0.5 * (eta_jk_old / tuning$E_rev[5])^2 - 
                                        log(tuning$E_rev[6]) - 0.5 * (theta_jk_old / tuning$E_rev[6])^2) -
                  (U_kj_new == 1) * (-log(tuning$E_rev[5]) - 0.5 * (eta_kj_new / tuning$E_rev[5])^2 - 
                                        log(tuning$E_rev[6]) - 0.5 * (theta_kj_new / tuning$E_rev[6])^2)
               ratio_MH    = exp((ratio_likel + ratio_prior) / temp + ratio_prop)
               
               # accept proposal with probability of min(1, MH ratio)
               if (is.nan(ratio_MH)) ratio_MH = 0
               if (runif(1) < min(1, ratio_MH))
               {
                  A0[j, k] = A0_jk_new
                  A0[k, j] = A0_kj_new
                  A1[j, k] = A1_jk_new
                  A1[k, j] = A1_kj_new
                  U[j, k]  = U_jk_new
                  U[k, j]  = U_kj_new
                  eta[j, k]    = eta_jk_new
                  eta[k, j]    = eta_kj_new
                  theta[j, k]  = theta_jk_new
                  theta[k, j]  = theta_kj_new
                  alpha0[j, k] = alpha0_jk_new
                  alpha0[k, j] = alpha0_kj_new
                  alpha1[j, k] = alpha1_jk_new
                  alpha1[k, j] = alpha1_kj_new
                  beta0[j, k]  = beta0_jk_new 
                  beta0[k, j]  = beta0_kj_new 
                  beta1[j, k]  = beta1_jk_new
                  beta1[k, j]  = beta1_kj_new 
                  delta0[j]    = delta0_j_new 
                  delta0[k]    = delta0_k_new 
                  delta1[j]    = delta1_j_new 
                  delta1[k]    = delta1_k_new 
                  gamma0[j]    = gamma0_j_new 
                  gamma0[k]    = gamma0_k_new 
                  gamma1[j]    = gamma1_j_new 
                  gamma1[k]    = gamma1_k_new
                  logitPi0[ , j] = logitPi0_j_new
                  logitPi0[ , k] = logitPi0_k_new
                  logitPi1[ , j] = logitPi1_j_new
                  logitPi1[ , k] = logitPi1_k_new
                  logMu0[ , j]   = logMu0_j_new
                  logMu0[ , k]   = logMu0_k_new
                  logMu1[ , j]   = logMu1_j_new
                  logMu1[ , k]   = logMu1_k_new
               }
            }
         }
      }
   }
   
   # update tau
   tau[1] = rgamma(1, shape = (0.5 * (sum(A0) + sum(A1 * (1 - A0))) + priors$tau[1] + temp - 1) / temp,
                   rate = (0.5 * (sum(alpha0 * alpha0) + sum(alpha1 * alpha1 * A1 * (1 - A0))) + priors$tau[2]) / temp)
   tau[2] = rgamma(1, shape = (0.5 * (sum(A0) + sum(A1 * (1 - A0))) + priors$tau[1] + temp - 1) / temp, 
                   rate = (0.5 * (sum(beta0 * beta0) + sum(beta1 * beta1 * A1 * (1 - A0))) + priors$tau[2]) / temp)
   tau[3] = rgamma(1, shape = (p + priors$tau[1] + temp - 1) / temp, 
                   rate = (0.5 * sum(delta0 * delta0 + delta1 * delta1) + priors$tau[2]) / temp)
   tau[4] = rgamma(1, shape = (p + priors$tau[1] + temp - 1) / temp, 
                   rate = (0.5 * sum(gamma0 * gamma0 + gamma1 * gamma1) + priors$tau[2]) / temp)
   tau[5] = rgamma(1, shape = (0.5 * sum(U, na.rm = TRUE) + priors$tau[1] + temp - 1) / temp, 
                   rate = (0.5 * sum(eta * eta, na.rm = TRUE) + priors$tau[2]) / temp)
   tau[6] = rgamma(1, shape = 0.5 * sum(U, na.rm = TRUE) + priors$tau[1] + temp - 1, 
                   rate = (0.5 * sum(theta * theta, na.rm = TRUE) + priors$tau[2]) / temp)
   
   # update rho
   rho[1] = rbeta(1, shape1 = (sum(U, na.rm = TRUE) + priors$rho[1] + temp - 1) / temp, 
                  shape2 = (sum(A0 * A1) - sum(U, na.rm = TRUE) + priors$rho[2] + temp - 1) / temp)
   rho[2] = rbeta(1, shape1 = (sum(D) + priors$rho[1] + temp - 1) / temp, 
                  shape2 = (p * (p - 1) - sum(D) + priors$rho[2] + temp - 1) / temp)
   rho[3] = rbeta(1, shape1 = (sum(A0) + priors$rho[1] + temp - 1) / temp, 
                  shape2 = (p * (p - 1) - sum(A0) + priors$rho[2] + temp - 1) / temp)
   
   # return updated parameters and acceptance record
   result = list()
   result$param = list()
   result$param$alpha1   = alpha1
   result$param$alpha0   = alpha0
   result$param$beta1    = beta1
   result$param$beta0    = beta0
   result$param$delta1   = delta1
   result$param$delta0   = delta0
   result$param$gamma1   = gamma1
   result$param$gamma0   = gamma0
   result$param$psi1     = psi1
   result$param$psi0     = psi0
   result$param$zeta     = eta            # change label
   result$param$eta      = theta          # change label
   result$param$U        = U
   result$param$E1       = A1             # change label
   result$param$D        = D
   result$param$E0       = A0             # change label
   result$param$tau      = tau
   result$param$rho      = rho
   result$param$logitPi1 = logitPi1
   result$param$logitPi0 = logitPi0
   result$param$logMu1   = logMu1
   result$param$logMu0   = logMu0
   result$acceptance = acceptance
   
   return(result)
}

# compute log of joint density of differential network models
log_dDN = function(x0, x1, param, priors, temp = 1)
{
   # the number of nodes
   p = ncol(x1)
   
   # given parameters
   alpha1   = param$alpha1
   alpha0   = param$alpha0
   beta1    = param$beta1
   beta0    = param$beta0
   delta1   = param$delta1
   delta0   = param$delta0
   gamma1   = param$gamma1
   gamma0   = param$gamma0
   psi1     = param$psi1
   psi0     = param$psi0
   eta      = param$zeta         # change label
   theta    = param$eta          # change label
   U        = param$U
   A1       = param$E1           # change label
   D        = param$D
   A0       = param$E0           # change label
   tau      = param$tau
   rho      = param$rho
   logitPi1 = param$logitPi1
   logitPi0 = param$logitPi0
   logMu1   = param$logMu1
   logMu0   = param$logMu0
   
   # calculate the likelihood part
   llik= 0
   for (j in 1 : p)
   {
      llik = llik + sum(llik_ZINBBN_j(x1[ , j], logitPi1[ , j], logMu1[ , j], psi1[j])) + 
         sum(llik_ZINBBN_j(x0[ , j], logitPi0[ , j], logMu0[ , j], psi0[j]))
   }
   
   # calculate the prior part
   lprior = sum(dnorm(alpha1[A1 * (1 - A0) == 1], mean = 0, sd = sqrt(1 / tau[1]), log = TRUE)) + 
      sum(dnorm(alpha0[A0 == 1], mean = 0, sd = sqrt(1 / tau[1]), log = TRUE)) + 
      sum(dnorm(beta1[A1 * (1 - A0) == 1], mean = 0, sd = sqrt(1 / tau[2]), log = TRUE)) + 
      sum(dnorm(beta0[A0 == 1], mean = 0, sd = sqrt(1 / tau[2]), log = TRUE)) + 
      sum(dnorm(delta1, mean = 0, sd = sqrt(1 / tau[3]), log = TRUE)) + 
      sum(dnorm(delta0, mean = 0, sd = sqrt(1 / tau[3]), log = TRUE)) +
      sum(dnorm(gamma1, mean = 0, sd = sqrt(1 / tau[4]), log = TRUE)) + 
      sum(dnorm(gamma0, mean = 0, sd = sqrt(1 / tau[4]), log = TRUE)) + 
      sum(dgamma(psi1, shape = priors$psi[1], rate = priors$psi[2], log = TRUE)) + 
      sum(dgamma(psi0, shape = priors$psi[1], rate = priors$psi[2], log = TRUE)) + 
      sum(dnorm(na.omit(eta[U == 1]), mean = 0, sd = sqrt(1 / tau[5]), log = TRUE)) + 
      sum(dnorm(na.omit(theta[U == 1]), mean = 0, sd = sqrt(1 / tau[6]), log = TRUE)) + 
      sum(U == 1, na.rm = TRUE) * log(rho[1]) + sum(U == 0, na.rm = TRUE) * log(1 - rho[1]) + 
      sum(D == 1) * log(rho[2]) + (sum(D == 0) - p) * log(1 - rho[2]) +
      sum(A0 == 1) * log(rho[3]) + (sum(A0 == 0) - p) * log(1 - rho[3]) +
      dgamma(tau[1], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) +
      dgamma(tau[2], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) +
      dgamma(tau[3], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) +
      dgamma(tau[4], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) +
      dgamma(tau[5], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) +
      dgamma(tau[6], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) +
      dbeta(rho[1], shape1 = priors$rho[1], shape2 = priors$rho[2], log = TRUE) + 
      dbeta(rho[2], shape1 = priors$rho[1], shape2 = priors$rho[2], log = TRUE) + 
      dbeta(rho[3], shape1 = priors$rho[1], shape2 = priors$rho[2], log = TRUE) 

   # return log-joint density
   return((llik + lprior) / temp)
}

# evaluate log-likelihood of each observation for node j of ZINBBN model
llik_ZINBBN_j = function(x_j, logitPi_j, logMu_j, psi_j)
{
   # calculate pi and mu for node j
   pi_j = exp(logitPi_j) / (1 + exp(logitPi_j))
   mu_j = exp(logMu_j)
   pi_j[is.nan(pi_j)] = 1
   mu_j[mu_j == Inf]  = .Machine$double.xmax
   
   # evaluate and return log-likelihood of each observation for node j
   llik_j = rep(0, length(x_j))
   llik_j[x_j == 0] = log(pi_j[x_j == 0] + (1 - pi_j[x_j == 0]) * dnbinom(0, size = 1 / psi_j, mu = mu_j[x_j == 0], log = FALSE))
   llik_j[x_j > 0]  = log(1 - pi_j[x_j > 0]) + dnbinom(x_j[x_j > 0], size = 1 / psi_j, mu = mu_j[x_j > 0], log = TRUE)
   return(llik_j)
}