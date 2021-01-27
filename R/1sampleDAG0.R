#' Generate data from one-sample DAG0 with given DAG and model parameters
#'
#' @param n the sample size.
#' @param E an adjacency matrix E of a DAG such that e_\{jk\} = 1 if k -> j. 
#' @param alpha a matrix for the parameter alpha of one-sample DAG0.
#' @param beta a matrix for the parameter beta of one-sample DAG0. 
#' @param delta a vector for the parameter delta of one-sample DAG0. 
#' @param gamma a vector for the parameter gamma of one-sample DAG0. 
#' @param psi a vector for the parameter psi of one-sample DAG0. 
#'
#' @return a matrix containing data generated from one-sample DAG0 with given DAG and model parameters.
#' @export
#'
#' @examples
#' library(BayesDAG0)
#' library(igraph)
#' 
#' # set random seed
#' set.seed(1)
#' 
#' # generate a random DAG E
#' p   = 10       # the number of variables 
#' n_e = p        # the number of edges
#' E_true = matrix(0, p, p)
#' while (sum(E_true == 1) < n_e)
#' {
#'    id_edge = matrix(sample(1 : p, 2), ncol = 2)
#'    E_true[id_edge] = 1
#'    g_true = graph_from_adjacency_matrix(t(E_true))
#'    if (!(is_dag(g_true)))
#'       E_true[id_edge] = 0
#' }
#' 
#' # generate model parameters for one-sample DAG0
#' alpha_true = matrix(0, p, p)
#' alpha_true[E_true == 1] = runif(n_e, 0.3, 1)
#' beta_true  = matrix(0, p, p)
#' beta_true[E_true == 1]  = runif(n_e, -1, -0.3)
#' delta_true = runif(p, -2, -1)
#' gamma_true = runif(p, 1, 2)
#' psi_true   = 1 / runif(p, 1, 5)
#' 
#' # generate data from one-sample DAG0 with given DAG and model parameters
#' dat = generate_data_DAG0(n = 100, E_true, alpha_true, beta_true, delta_true, gamma_true, psi_true)
generate_data_DAG0 = function(n, E, alpha, beta, delta, gamma, psi)
{
   # generate data from one-sample DAG0
   dat = matrix(0, n, p)
   G   = graph_from_adjacency_matrix(t(E))
   order_nodes = as_ids(topo_sort(G))
   for (j in order_nodes)
   {
      pi = exp(dat %*% alpha[j, ] + delta[j])
      pi = pi / (1 + pi)
      pi[is.nan(pi)] = 1
      mu = exp(dat %*% beta[j, ] + gamma[j])
      dat[ , j] = rnbinom(n, size = 1 / psi[j], mu = mu) * (1 - rbinom(n, size = 1, prob = pi))
   }
   
   return(dat)
}

#' Implementation of the parallel-tempered MCMC for one-sample DAG0
#'
#' @param data a matrix containing data.
#' @param starting a list with each tag corresponding to a parameter name. 
#' Valid tags are 'E', 'alpha', 'beta', 'delta', 'gamma', 'psi', 'tau', and 'rho'. 
#' The value portion of each tag is the parameters' starting values for MCMC.
#' @param tuning a list with each tag corresponding to a parameter name. 
#' Valid tags are 'E', 'E_rev', 'alpha', 'beta', 'delta', 'gamma', and 'psi'. 
#' The value portion of each tag defines the standard deviations of Normal proposal distributions for Metropolis sampler for each parameter. 
#' The tag 'E' corresponds to birth or death of an edge, while 'E_rev' indicates reversal of an edge.  
#' The value portion of both tags should consist of the standard deviations of Gaussian proposals for alpha, beta, delta, and gamma for an update of E.
#' @param priors a list with each tag corresponding to a parameter name. 
#' Valid tags are 'psi', 'tau', and 'rho'. The value portion of each tag defines the hyperparameters for one-sample DAG0.
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
#' The value portion of the tag 'samples' gives MCMC samples from the posterior distribution of the parameters for one-sample DAG0.
#' The value portion of the tag 'acceptance' shows the Metropolis sampling acceptance percents for alpha, beta, delta, gamma, and psi. 
#' @export
#'
#' @examples
#' library(BayesDAG0)
#' library(igraph)
#' library(pscl)
#' 
#' # set random seed
#' set.seed(1)
#' 
#' # generate a random DAG E
#' p   = 10       # the number of variables 
#' n_e = p        # the number of edges
#' E_true = matrix(0, p, p)
#' while (sum(E_true == 1) < n_e)
#' {
#'    id_edge = matrix(sample(1 : p, 2), ncol = 2)
#'    E_true[id_edge] = 1
#'    g_true = graph_from_adjacency_matrix(t(E_true))
#'    if (!(is_dag(g_true)))
#'       E_true[id_edge] = 0
#' }
#' 
#' # generate model parameters for one-sample DAG0
#' alpha_true = matrix(0, p, p)
#' alpha_true[E_true == 1] = runif(n_e, 0.3, 1)
#' beta_true  = matrix(0, p, p)
#' beta_true[E_true == 1]  = runif(n_e, -1, -0.3)
#' delta_true = runif(p, -2, -1)
#' gamma_true = runif(p, 1, 2)
#' psi_true   = 1 / runif(p, 1, 5)
#' 
#' # generate data from one-sample DAG0 with given DAG and model parameters
#' dat = generate_data_DAG0(n = 100, E_true, alpha_true, beta_true, delta_true, gamma_true, psi_true)
#' 
#' ### conduct posterior inference for one-sample DAG0
#' # starting values for MCMC 
#' starting = list()
#' starting$E     = matrix(0, p, p)
#' starting$alpha = matrix(0, p, p)
#' starting$beta  = matrix(0, p, p)
#' starting$delta = rep(0, p)
#' starting$gamma = rep(0, p)
#' starting$psi   = rep(0, p)
#' for (j in 1 : p)
#' {
#'    mle = zeroinfl(dat[ , j] ~ 1 | 1, dist = "negbin")
#'    starting$delta[j] = mle$coefficients$zero
#'    starting$gamma[j] = mle$coefficients$count
#'    starting$psi[j]   = 1 / mle$theta
#' }
#' 
#' starting$tau = c(10, 10, 1, 1)
#' starting$rho = 1 / p
#' 
#' # sd's of the Metropolis sampler Normal proposal distribution
#' tuning = list()
#' tuning$E     = c(0.05, 0.05, 1.2, 0.3)
#' tuning$E_rev = c(0.5, 0.5, 1.2, 0.3) 
#' tuning$alpha = 1.2  
#' tuning$beta  = 0.4  
#' tuning$delta = 1.2  
#' tuning$gamma = 0.3 
#' tuning$psi   = 0.8  
#' 
#' # hyperparameters for one-sample DAG0
#' priors = list()
#' priors$psi = c(1, 1)
#' priors$tau = c(1, 1)
#' priors$rho = c(1, 1)
#' 
#' # run the parallel-tempered MCMC for one-sample DAG0
#' out = mcmc_1sampleDAG0(dat, starting, tuning, priors, n_sample = 2000, n_burnin = 1000, do_par = TRUE)
#' 
#' # report Metropolis sampling acceptance rates
#' out$acceptance
#' 
#' # graph estimate based on median probability model
#' E_est  = (apply(out$samples$E, c(1, 2), mean) > 0.5) + 0
#' 
#' # report contingency table for edge selection
#' table(E_est, E_true)
mcmc_1sampleDAG0 = function(data, starting, tuning, priors, n_sample = 5000, n_burnin = 2500, 
                            n_chain = 10, prob_swap = 0.1, temperature = NULL, do_par = FALSE, n_core = n_chain,
                            verbose = TRUE, n_report = 100)
{
   # sample size and dimension
   n = nrow(data)
   p = ncol(data)
   
   # default temperatures for the parallel tempering 
   if (is.null(temperature))
      temperature = exp(seq(0, log(1.5), length.out = n_chain))
   
   # sd's of the Metropolis sampler Normal proposal distribution for each chain
   list_tuning = list()
   for (m in 1 : n_chain)
   {
      list_tuning[[m]] = list()
      list_tuning[[m]]$E     = tuning$E
      list_tuning[[m]]$E_rev = tuning$E_rev
      list_tuning[[m]]$alpha = tuning$alpha * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$beta  = tuning$beta  * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$delta = tuning$delta * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$gamma = tuning$gamma * temperature[n_chain]^((m - 1) / (n_chain - 1))
      list_tuning[[m]]$psi   = tuning$psi   * temperature[n_chain]^((m - 1) / (n_chain - 1))
   }

   # calculate logitPi and logLambda with starting values
   starting$logitPi = tcrossprod(data, starting$alpha) + matrix(starting$delta, n, p, byrow = TRUE)
   starting$logMu   = tcrossprod(data, starting$beta) + matrix(starting$gamma, n, p, byrow = TRUE)

   # initialize parameters 
   param = list()
   for (m in 1 : n_chain)
   {
      param[[m]] = starting
   }
   
   # initialize MCMC samples for the cold chain
   mcmc_samples = list()
   mcmc_samples$E     = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$alpha = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$beta  = array(NA, dim = c(p, p, n_sample))
   mcmc_samples$delta = matrix(NA, p, n_sample)
   mcmc_samples$gamma = matrix(NA, p, n_sample)
   mcmc_samples$psi   = matrix(NA, p, n_sample)
   mcmc_samples$tau   = matrix(NA, 4, n_sample)
   mcmc_samples$rho   = rep(NA, n_sample)
   
   # initialize acceptance indicators for the cold chain
   acceptance = list()
   acceptance$alpha = array(NA, dim = c(p, p, n_sample))
   acceptance$beta  = array(NA, dim = c(p, p, n_sample))
   acceptance$delta = matrix(NA, p, n_sample)
   acceptance$gamma = matrix(NA, p, n_sample)
   acceptance$psi   = matrix(NA, p, n_sample)
   
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
               gibbs_1sampleDAG0(param[[m]], list_tuning[[m]], priors, data, temperature[m])
            }
         } else               # update chains sequentially
         {
            out = foreach(m = 1 : n_chain, .packages = "igraph") %do% {
               gibbs_1sampleDAG0(param[[m]], list_tuning[[m]], priors, data, temperature[m])
            }
         }
         
         for (m in 1 : n_chain)
         {
            param[[m]] = out[[m]]$param
            
            if (m == 1)
            {
               # save MCMC samples for the cold chain
               mcmc_samples$E[ , , t]     = param[[1]]$E
               mcmc_samples$alpha[ , , t] = param[[1]]$alpha
               mcmc_samples$beta[ , , t]  = param[[1]]$beta
               mcmc_samples$delta[ , t]   = param[[1]]$delta
               mcmc_samples$gamma[ , t]   = param[[1]]$gamma
               mcmc_samples$psi[ , t]     = param[[1]]$psi
               mcmc_samples$tau[ , t]     = param[[1]]$tau
               mcmc_samples$rho[t]        = param[[1]]$rho
               
               # save acceptance rates for the cold chain
               acceptance$alpha[ , , t] = out[[1]]$acceptance$alpha
               acceptance$beta[ , , t]  = out[[1]]$acceptance$beta
               acceptance$delta[ , t]   = out[[1]]$acceptance$delta
               acceptance$gamma[ , t]   = out[[1]]$acceptance$gamma
               acceptance$psi[ , t]     = out[[1]]$acceptance$psi
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
         
         # calculate MH ratio for the swapping
         ratio_swap = log_dZINBBN(data, param[[m1]], priors, temperature[m2]) +
            log_dZINBBN(data, param[[m2]], priors, temperature[m1]) -
            log_dZINBBN(data, param[[m1]], priors, temperature[m1]) -
            log_dZINBBN(data, param[[m2]], priors, temperature[m2])
         ratio_swap = exp(ratio_swap)
         
         # accept the swapping
         if (is.nan(ratio_swap)) ratio_swap = 0
         if (runif(1) < min(1, ratio_swap))
         {
            param_m1    = param[[m1]]
            param[[m1]] = param[[m2]]
            param[[m2]] = param_m1
         }
         
         # save MCMC samples for the  cold chain
         mcmc_samples$E[ , , t]     = param[[1]]$E
         mcmc_samples$alpha[ , , t] = param[[1]]$alpha
         mcmc_samples$beta[ , , t]  = param[[1]]$beta
         mcmc_samples$delta[ , t]   = param[[1]]$delta
         mcmc_samples$gamma[ , t]   = param[[1]]$gamma
         mcmc_samples$psi[ , t]     = param[[1]]$psi
         mcmc_samples$tau[ , t]     = param[[1]]$tau
         mcmc_samples$rho[t]        = param[[1]]$rho
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
   # 2. Metropolis sampler acceptance rates for alpha, beta, delta, gamma, and psi
   results = list()
   results$samples = list(E     = mcmc_samples$E[ , , (n_burnin + 1) : n_sample],
                          alpha = mcmc_samples$alpha[ , , (n_burnin + 1) : n_sample],
                          beta  = mcmc_samples$beta[ , , (n_burnin + 1) : n_sample],
                          delta = mcmc_samples$delta[ , (n_burnin + 1) : n_sample],
                          gamma = mcmc_samples$gamma[ , (n_burnin + 1) : n_sample],
                          psi   = mcmc_samples$psi[ , (n_burnin + 1) : n_sample],
                          tau   = mcmc_samples$tau[ , (n_burnin + 1) : n_sample],
                          rho   = mcmc_samples$rho[(n_burnin + 1) : n_sample])
   results$acceptance = list(alpha = 100 * mean(acceptance$alpha, na.rm = TRUE),
                             beta  = 100 * mean(acceptance$beta, na.rm = TRUE),
                             delta = 100 * mean(acceptance$delta, na.rm = TRUE),
                             gamma = 100 * mean(acceptance$gamma, na.rm = TRUE),
                             psi   = 100 * mean(acceptance$psi, na.rm = TRUE))
   class(results) = "1sampleDAG0"
   
   return(results)
}

#' Implementation of Gibbs step for one-sample DAG0
#'
#' @param param a list with each tag corresponding to a parameter name. 
#' Valid tags are 'E', 'alpha', 'beta', 'delta', 'gamma', 'psi', 'tau', and 'rho'. 
#' The value portion of each tag is the parameter values to be updated.
#' @param tuning a list with each tag corresponding to a parameter name. 
#' Valid tags are 'E', 'E_rev', 'alpha', 'beta', 'delta', 'gamma', and 'psi'. 
#' The value portion of each tag defines the standard deviations of Normal proposal distributions for Metropolis sampler for each parameter. 
#' The tag 'E' corresponds to birth or death of an edge, while 'E_rev' indicates reversal of an edge.  
#' The value portion of both tags should consist of the standard deviations of Gaussian proposals for alpha, beta, delta, and gamma for an update of E.
#' @param priors a list with each tag corresponding to a parameter name. 
#' Valid tags are 'psi', 'tau', and 'rho'. The value portion of each tag defines the hyperparameters for one-sample DAG0.
#' @param data a matrix containing data.
#' @param temp a temperature.
#' @param prob_rev the probability of edge-reversal.
#'
#' @return a list with the tags 'param' and 'acceptance'. 
#' The value portion of the tag 'param' gives parameter values after perfomring Gibbs step for one-sample DAG0.
#' The value portion of the tag 'acceptance' shows the acceptance indicators for alpha, beta, delta, gamma, and psi.
#' @export
#' 
#' @examples
#' 
gibbs_1sampleDAG0 = function(param, tuning, priors, data, temp = 1, prob_rev = 0.5)
{
   # parameters to be updated
   alpha   = param$alpha
   beta    = param$beta
   delta   = param$delta
   gamma   = param$gamma
   psi     = param$psi
   A       = param$E          # change label
   tau     = param$tau
   rho     = param$rho
   logitPi = param$logitPi
   logMu   = param$logMu

   # update alpha 
   update1 = update_alpha(data, alpha, psi, A, tau, logitPi, logMu, tuning$alpha^(-2), temp)
   alpha   = update1$alpha
   logitPi = update1$logitPi
   
   # update beta
   update2 = update_beta(data, beta, psi, A, tau, logitPi, logMu, tuning$beta^(-2), temp)
   beta    = update2$beta
   logMu   = update2$logMu
   
   # update delta
   update3 = update_delta(data, delta, psi, tau, logitPi, logMu, tuning$delta^(-2), temp)
   delta   = update3$delta
   logitPi = update3$logitPi
   
   # update gamma 
   update4 = update_gamma(data, gamma, psi, tau, logitPi, logMu, tuning$gamma^(-2), temp)
   gamma   = update4$gamma
   logMu   = update4$logMu
   
   # update psi
   update5 = update_psi(data, psi, logitPi, logMu, tuning$psi^(-2), priors, temp)
   psi     = update5$psi
   
   # update A
   if (runif(1) > prob_rev)
      update6 = update_A_each(data, alpha, beta, delta, gamma, psi, A, tau, rho, logitPi, logMu, tuning$E^(-2), temp)
   else
      update6 = update_A_rev(data, alpha, beta, delta, gamma, psi, A, tau, logitPi, logMu, tuning$E_rev^(-2), temp)
   A       = update6$A
   alpha   = update6$alpha
   beta    = update6$beta
   delta   = update6$delta
   gamma   = update6$gamma
   logitPi = update6$logitPi
   logMu   = update6$logMu
   
   # update tau
   tau = update_tau(alpha, beta, delta, gamma, A, priors, temp)
   
   # uphdate rho
   rho = update_rho(A, priors, temp)
   
   # return the updated parameters and the acceptance record
   result = list()
   result$param      = list()
   result$acceptance = list()
   result$param$alpha   = alpha
   result$param$beta    = beta
   result$param$delta   = delta
   result$param$gamma   = gamma
   result$param$psi     = psi
   result$param$E       = A         # change label
   result$param$tau     = tau
   result$param$rho     = rho
   result$param$logitPi = logitPi
   result$param$logMu   = logMu
   result$acceptance$alpha = update1$accept
   result$acceptance$beta  = update2$accept
   result$acceptance$delta = update3$accept
   result$acceptance$gamma = update4$accept
   result$acceptance$psi   = update5$accept
   
   return(result)
}

# compute log-joint density of ZINBBN models
log_dZINBBN = function(data, param, priors, temp = 1)
{
   p = ncol(data)
   
   # specify parameters 
   alpha   = param$alpha
   beta    = param$beta
   delta   = param$delta
   gamma   = param$gamma
   psi     = param$psi
   A       = param$E
   tau     = param$tau
   rho     = param$rho
   logitPi = param$logitPi
   logMu   = param$logMu
   
   # calculate the likelihood
   llik = 0
   for (j in 1 : p)
      llik = llik + sum(llik_ZINBBN_j(data[ , j], logitPi[ , j], logMu[ , j], psi[j]))
   
   # calculate the prior density
   lprior = sum(dnorm(alpha[A == 1], mean = 0, sd = sqrt(1 / tau[1]), log = TRUE)) + 
      sum(dnorm(beta[A == 1], mean = 0, sd = sqrt(1 / tau[2]), log = TRUE)) + 
      sum(dnorm(delta, mean = 0, sd = sqrt(1 / tau[3]), log = TRUE)) + 
      sum(dnorm(gamma, mean = 0, sd = sqrt(1 / tau[4]), log = TRUE)) + 
      sum(dgamma(psi, shape = priors$psi[1], rate = priors$psi[2], log = TRUE)) + 
      sum(A) * log(rho) + (p * (p - 1) - sum(A)) * log(1 - rho) + 
      dgamma(tau[1], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) + 
      dgamma(tau[2], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) + 
      dgamma(tau[3], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) + 
      dgamma(tau[4], shape = priors$tau[1], rate = priors$tau[2], log = TRUE) + 
      dbeta(rho, shape1 = priors$rho[1], shape2 = priors$rho[2], log = TRUE)
   
   return((llik + lprior) / temp)
}

# update each element of alpha through Metropolis-Hastings step
update_alpha = function(x, alpha, psi, A, tau, logitPi, logMu, phi_alpha, temp = 1)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = matrix(NA, p, p)
   
   for (j in 1 : p)
   {
      # get data, logit(pi), log(mu), and psi for node j
      x_j = x[ , j]
      logitPi_j = logitPi[ , j]
      logMu_j   = logMu[ , j]
      psi_j     = psi[j]
      
      for (k in 1 : p)
      {
         # if j = k, do not sample 
         if (j == k) next
         
         if (A[j, k] == 1)
         {
            accept[j, k] = 0
            
            # current value of \alpha_{jk}
            alpha_old = alpha[j, k]
            
            # propose new value of \alpha_{jk} using Normal proposal distribution
            alpha_new   = rnorm(1, mean = alpha_old, sd = sqrt(1 / phi_alpha))
            logitPi_new = logitPi_j + x[ , k] * (alpha_new - alpha_old)
            
            # calculate MH ratio
            llik_old = llik_ZINBBN_j(x_j, logitPi_j, logMu_j, psi_j)
            llik_new = llik_ZINBBN_j(x_j, logitPi_new, logMu_j, psi_j)
            ratio_MH = exp((sum(llik_new - llik_old) - 0.5 * tau[1] * (alpha_new * alpha_new - alpha_old * alpha_old)) / temp)
         
            # accept the proposed value with probability of min(1, MH ratio)
            if (is.nan(ratio_MH)) ratio_MH = 0
            if (runif(1) < min(1, ratio_MH))
            {
               alpha[j, k]  = alpha_new
               logitPi_j    = logitPi_new   # update logit(pi) for node j 
               accept[j, k] = 1             # 1 if proposal accepted
            }
         }
      }
      
      # save the updated logit(pi) for the node j
      logitPi[ , j] = logitPi_j
   }
   
   # return the updated alpha and logit(pi), and the acceptance indicators 
   return(list(alpha   = alpha, 
               logitPi = logitPi,
               accept  = accept))
}

# update each element of beta through Metropolis-Hastings step
update_beta = function(x, beta, psi, A, tau, logitPi, logMu, phi_beta, temp = 1)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = matrix(NA, p, p)
   
   for (j in 1 : p)
   {
      # get data, logit(pi), log(mu), and psi for node j
      x_j = x[ , j]
      logitPi_j = logitPi[ , j]
      logMu_j   = logMu[ , j]
      psi_j     = psi[j]
      
      for (k in 1 : p)
      {
         # if j = k, do not sample 
         if (j == k) next
         
         if(A[j, k] == 1)
         {
            accept[j, k] = 0
            
            # current value of \beta_{jk}
            beta_old = beta[j, k]
      
            # propose new value of \beta_{jk} using Normal proposal distribution
            beta_new  = rnorm(1, mean = beta_old, sd = sqrt(1 / phi_beta))
            logMu_new = logMu_j + x[ , k] * (beta_new - beta_old)
            
            # calculate MH ratio
            llik_old = llik_ZINBBN_j(x_j, logitPi_j, logMu_j, psi_j)
            llik_new = llik_ZINBBN_j(x_j, logitPi_j, logMu_new, psi_j)
            ratio_MH = exp((sum(llik_new - llik_old) - 0.5 * tau[2] * (beta_new * beta_new - beta_old * beta_old)) / temp)
            
            # accept the proposed value with probability of min(1, MH ratio)
            if (is.nan(ratio_MH)) ratio_MH = 0
            if (runif(1) < min(1, ratio_MH))
            {
               beta[j, k]   = beta_new
               logMu_j      = logMu_new   # update log(mu) for node j
               accept[j, k] = 1           # 1 if proposal accepted
            }
         }
      }
      
      # save the updated log(mu) for the node j
      logMu[ , j] = logMu_j
   }
   
   # return the updated beta and log(mu), and the acceptance indicators 
   return(list(beta      = beta,
               logMu = logMu,
               accept    = accept))
}

# update each element of delta through Metropolis-Hastings step
update_delta = function(x, delta, psi, tau, logitPi, logMu, phi_delta, temp = 1)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = rep(0, p)
   
   for (j in 1 : p)
   {
      # get data, logit(pi), log(mu), and psi for node j
      x_j       = x[ , j]
      logitPi_j = logitPi[ , j]
      logMu_j   = logMu[ , j]
      psi_j     = psi[j]
      
      # current value of \delta_j
      delta_old   = delta[j]
      
      # propose new value of \delta_j using Normal proposal distribution
      delta_new   = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_delta))
      logitPi_new = logitPi_j + (delta_new - delta_old)
      
      # calculate MH ratio
      llik_old = llik_ZINBBN_j(x_j, logitPi_j, logMu_j, psi_j)
      llik_new = llik_ZINBBN_j(x_j, logitPi_new, logMu_j, psi_j)
      ratio_MH = exp((sum(llik_new - llik_old) - 0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old)) / temp)
      
      # accept the proposed value with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
      if (runif(1) < min(1, ratio_MH))
      {
         delta[j]      = delta_new
         logitPi[ , j] = logitPi_new   # update and save logit(pi) for node j 
         accept[j]     = 1             # 1 if proposal accepted
      }
   }
   
   # return the updated delta and logit(pi), and the acceptance indicators
   return(list(delta   = delta,
               logitPi = logitPi,
               accept  = accept))
}

# update each element of gamma through Metropolis-Hastings step
update_gamma = function(x, gamma, psi, tau, logitPi, logMu, phi_gamma, temp = 1)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = rep(0, p)
   
   for (j in 1 : p)
   {
      # get data, logit(pi), log(mu), and psi for node j
      x_j       = x[ , j]
      logitPi_j = logitPi[ , j]
      logMu_j   = logMu[ , j]
      psi_j     = psi[j]
      
      # current value of \gamma_j
      gamma_old = gamma[j]
      
      # propose new value of \gamma_j using Normal proposal distribution
      gamma_new     = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_gamma))
      logMu_new = logMu_j + (gamma_new - gamma_old)
      
      # calculate MH ratio
      llik_old = llik_ZINBBN_j(x_j, logitPi_j, logMu_j, psi_j)
      llik_new = llik_ZINBBN_j(x_j, logitPi_j, logMu_new, psi_j)
      ratio_MH = exp((sum(llik_new - llik_old) - 0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old)) / temp)
      
      # accept the proposed value with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
      if (runif(1) < min(1, ratio_MH))
      {
         gamma[j]    = gamma_new
         logMu[ , j] = logMu_new   # update and save log(mu) for  node j
         accept[j]   = 1           # 1 if proposal accepted
      }
   }
   
   # return the updated gamma and log(mu), and the acceptance indicators
   return(list(gamma  = gamma,
               logMu  = logMu,
               accept = accept))
}

# update each element of psi through Metropolis-Hastings step
update_psi = function(x, psi, logitPi, logMu, phi_psi, priors, temp = 1)
{
   # get the number of nodes and initialize acceptance indicators
   p      = ncol(x)
   accept = rep(0, p)
   
   for (j in 1 : p)
   {
      # get data, logit(pi), and log(mu) for node j
      x_j       = x[ , j]
      logitPi_j = logitPi[ , j]
      logMu_j   = logMu[ , j]
      
      # current value of log(\psi_j)
      logPsi_old = log(psi[j])
      
      # propose new value of log(\psi_j) using Normal proposal distribution
      logPsi_new = rnorm(1, mean = logPsi_old, sd = sqrt(1 / phi_psi))
      
      # calculate MH ratio
      llik_old = llik_ZINBBN_j(x_j, logitPi_j, logMu_j, exp(logPsi_old))
      llik_new = llik_ZINBBN_j(x_j, logitPi_j, logMu_j, exp(logPsi_new))
      ratio_MH = exp((sum(llik_new - llik_old) + priors$psi[1] * (logPsi_new - logPsi_old) - priors$psi[2] * (exp(logPsi_new) - exp(logPsi_old))) / temp)
      
      # accept the proposed valu with probability of min(1, MH ratio)
      if (is.nan(ratio_MH)) ratio_MH = 0
      if (runif(1) < min(1, ratio_MH))
      {
         psi[j]    = exp(logPsi_new)   # transformed back to psi
         accept[j] = 1                 # 1 if proposal accepted
      }
   }
   
   # return the updated psi and acceptance indicators
   return(list(psi    = psi,
               accept = accept))
}


# update each element of A through Metropolis-Hastings step (propose addition or deletion of an edge)
# corresponding alpha, beta, delta, and gamma are jointly proposed with A 
update_A_each = function(x, alpha, beta, delta, gamma, psi, A, tau, rho, logitPi, logMu, phi_A, temp = 1)
{
   # get the number of nodes
   p = ncol(x)
   
   for (j in 1 : p)
   {
      # get data, logit(pi), log(mu), and psi for node j
      x_j       = x[ , j]
      logitPi_j = logitPi[ , j]
      logMu_j   = logMu[ , j]
      psi_j     = psi[j]
      
      for (k in 1 : p)
      {
         # if j = k, do not sample 
         if (j == k) next
         
         # current value of \alpha_{jk}, \beta_{jk}, \delta_j, and \gamma_j
         alpha_old = alpha[j, k]
         beta_old  = beta[j, k]
         delta_old = delta[j]
         gamma_old = gamma[j]
         
         # divide into two cases: 1. there is currently no edge k -> j, 2. there exist an edge k -> j now 
         if (A[j, k] == 0)
         {
            # if there is no edge k -> j, propose addition of the edge unless it makes a cycle
            A[j, k] = A_new = 1
            graph   = graph_from_adjacency_matrix(A)
            A[j, k] = 0
            if (!is_dag(graph)) next
            
            # propose new value of \alpha_{jk}, \beta_{jk}, \delta_j, and \gamma_j 
            alpha_new = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
            beta_new  = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[2]))
            delta_new = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_A[3]))
            gamma_new = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_A[4]))
            logitPi_new = logitPi_j + x[ , k] * (alpha_new - alpha_old) + (delta_new - delta_old)
            logMu_new   = logMu_j + x[ , k] * (beta_new - beta_old) + (gamma_new - gamma_old)
            
            # calculate MH ratio
            llik_old   = llik_ZINBBN_j(x_j, logitPi_j, logMu_j, psi_j)
            llik_new   = llik_ZINBBN_j(x_j, logitPi_new, logMu_new, psi_j)
            ratio_post = sum(llik_new - llik_old) + 
               0.5 * log(tau[1]) - 0.5 * tau[1] * alpha_new * alpha_new + 
               0.5 * log(tau[2]) - 0.5 * tau[2] * beta_new * beta_new -
               0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old) - 
               0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old) +
               log(rho) - log(1 - rho)
            ratio_prop = -0.5 * log(phi_A[1]) + 0.5 * phi_A[1] * alpha_new * alpha_new - 
               0.5 * log(phi_A[2]) + 0.5 * phi_A[2] * beta_new * beta_new
            ratio_MH = exp(ratio_post / temp + ratio_prop)
         } 
         else
         {
            # if there is an edge k -> j, propose deletion of the edge
            A_new = 0
            
            # propose new value of \alpha_{jk}, \beta_{jk}, \delta_j, and \gamma_j 
            alpha_new = beta_new = 0
            delta_new = rnorm(1, mean = delta_old, sd = sqrt(1 / phi_A[3]))
            gamma_new = rnorm(1, mean = gamma_old, sd = sqrt(1 / phi_A[4]))
            logitPi_new = logitPi_j + x[ , k] * (alpha_new - alpha_old) + (delta_new - delta_old)
            logMu_new   = logMu_j + x[ , k] * (beta_new - beta_old) + (gamma_new - gamma_old)
            
            # calculate MH ratio
            llik_old   = llik_ZINBBN_j(x_j, logitPi_j, logMu_j, psi_j)
            llik_new   = llik_ZINBBN_j(x_j, logitPi_new, logMu_new, psi_j)
            ratio_post = sum(llik_new - llik_old) - 
               0.5 * log(tau[1]) + 0.5 * tau[1] * alpha_old * alpha_old - 
               0.5 * log(tau[2]) + 0.5 * tau[2] * beta_old * beta_old - 
               0.5 * tau[3] * (delta_new * delta_new - delta_old * delta_old) - 
               0.5 * tau[4] * (gamma_new * gamma_new - gamma_old * gamma_old) +
               log(1 - rho) - log(rho)
            ratio_prop = 0.5 * log(phi_A[1]) - 0.5 * phi_A[1] * alpha_old * alpha_old + 
               0.5 * log(phi_A[2]) - 0.5 * phi_A[2] * beta_old * beta_old
            ratio_MH = exp(ratio_post / temp + ratio_prop)
         }
         
         
         # accept the proposed values with probabiliof min(1, MH ratio) 
         if (is.nan(ratio_MH)) ratio_MH = 0
         if (runif(1) < min(1, ratio_MH))
         {
            A[j, k]     = A_new
            alpha[j, k] = alpha_new
            beta[j, k]  = beta_new
            delta[j]    = delta_new
            gamma[j]    = gamma_new
            logitPi_j   = logitPi_new   # update logit(pi) for node j
            logMu_j     = logMu_new     # update log(mu) for node j
         }
      }
      
      # save the updated logit(pi) and log(mu) for node j
      logitPi[ , j] = logitPi_j
      logMu[ , j]   = logMu_j
   }
   
   # return the updated A, alpha, beta, delta, gamma, logit(pi), and log(mu) 
   return(list(A       = A,
               alpha   = alpha,
               beta    = beta,
               delta   = delta,
               gamma   = gamma,
               logitPi = logitPi,
               logMu   = logMu))
}

# update A based on proposal of reversing an edge through Metropolis-Hastings step 
# corresponding alpha, beta, delta, and gamma are jointly proposed with A 
update_A_rev = function(x, alpha, beta, delta, gamma, psi, A, tau, logitPi, logMu, phi_A, temp = 1)
{
   # get indices of existing edges 
   id_edges = which(A == 1, arr.ind = TRUE)
   n_rev    = nrow(id_edges)
   
   # if there exists no edge, don't do anything
   if (n_rev > 0)
   {
      for (s in 1 : n_rev)
      {
         # index of an edge which will be reversed
         j = id_edges[s, 1]
         k = id_edges[s, 2]
         
         # propose reversal of the edge unless it makes a cycle
         A[j, k] = A_jk_new = 0
         A[k, j] = A_kj_new = 1
         graph   = graph_from_adjacency_matrix(A)
         A[j, k] = 1
         A[k, j] = 0
         if (!is_dag(graph)) next
         
         # get data, logit(pi), log(mu) and psi for reversal of the edge
         x_j       = x[ , j]
         x_k       = x[ , k]
         logitPi_j = logitPi[ , j]
         logitPi_k = logitPi[ , k]
         logMu_j   = logMu[ , j]
         logMu_k   = logMu[ , k]
         psi_j     = psi[j]
         psi_k     = psi[k]
         
         # current values of elements of alpha, beta, delta, and gamma corresponding to reversing
         alpha_jk_old = alpha[j, k]
         alpha_kj_old = alpha[k, j]
         beta_jk_old  = beta[j, k]
         beta_kj_old  = beta[k, j]
         delta_j_old  = delta[j]
         delta_k_old  = delta[k]
         gamma_j_old  = gamma[j]
         gamma_k_old  = gamma[k]
         
         # propose new values of elements of alpha, beta, delta, and gamma corresponding to reversing
         alpha_jk_new  = beta_jk_new  = 0
         alpha_kj_new  = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[1]))
         beta_kj_new   = rnorm(1, mean = 0, sd = sqrt(1 / phi_A[2]))
         delta_j_new   = rnorm(1, mean = delta_j_old, sd = sqrt(1 / phi_A[3]))
         delta_k_new   = rnorm(1, mean = delta_k_old, sd = sqrt(1 / phi_A[3]))
         gamma_j_new   = rnorm(1, mean = gamma_j_old, sd = sqrt(1 / phi_A[4]))
         gamma_k_new   = rnorm(1, mean = gamma_k_old, sd = sqrt(1 / phi_A[4]))
         logitPi_j_new = logitPi_j + x_k * (alpha_jk_new - alpha_jk_old) + (delta_j_new - delta_j_old)
         logitPi_k_new = logitPi_k + x_j * (alpha_kj_new - alpha_kj_old) + (delta_k_new - delta_k_old)
         logMu_j_new   = logMu_j + x_k * (beta_jk_new - beta_jk_old) + (gamma_j_new - gamma_j_old)
         logMu_k_new   = logMu_k + x_j * (beta_kj_new - beta_kj_old) + (gamma_k_new - gamma_k_old)
         
         # calculate MH ratio
         llik_j_old = llik_ZINBBN_j(x_j, logitPi_j, logMu_j, psi_j)
         llik_k_old = llik_ZINBBN_j(x_k, logitPi_k, logMu_k, psi_k)
         llik_j_new = llik_ZINBBN_j(x_j, logitPi_j_new, logMu_j_new, psi_j)
         llik_k_new = llik_ZINBBN_j(x_k, logitPi_k_new, logMu_k_new, psi_k)
         ratio_post = sum(llik_j_new - llik_j_old) + sum(llik_k_new - llik_k_old) - 
            0.5 * tau[1] * (alpha_kj_new * alpha_kj_new - alpha_jk_old * alpha_jk_old) - 
            0.5 * tau[2] * (beta_kj_new * beta_kj_new - beta_jk_old * beta_jk_old) -
            0.5 * tau[3] * (delta_j_new * delta_j_new - delta_j_old * delta_j_old + 
                               delta_k_new * delta_k_new - delta_k_old * delta_k_old) - 
            0.5 * tau[4] * (gamma_j_new * gamma_j_new - gamma_j_old * gamma_j_old +
                               gamma_k_new * gamma_k_new - gamma_k_old * gamma_k_old)
         ratio_prop = -0.5 * phi_A[1] * (alpha_jk_old * alpha_jk_old - alpha_kj_new * alpha_kj_new) - 
            0.5 * phi_A[2] * (beta_jk_old * beta_jk_old - beta_kj_new * beta_kj_new)
         ratio_MH   = exp(ratio_post / temp + ratio_prop)
         
         # accept the proposed value with probability of min(1, MH ratio)
         if (is.nan(ratio_MH)) ratio_MH = 0   
         if (runif(1) < min(1, ratio_MH))
         {
            A[j, k]     = A_jk_new
            A[k, j]     = A_kj_new
            alpha[j, k] = alpha_jk_new
            alpha[k, j] = alpha_kj_new
            beta[j, k]  = beta_jk_new
            beta[k, j]  = beta_kj_new
            delta[j]    = delta_j_new
            delta[k]    = delta_k_new
            gamma[j]    = gamma_j_new
            gamma[k]    = gamma_k_new
            logitPi[ , j] = logitPi_j_new
            logitPi[ , k] = logitPi_k_new   # update logit(pi)'s, following reversal of the edge
            logMu[ , j]   = logMu_j_new
            logMu[ , k]   = logMu_k_new     # update log(mu)'s, following reversal of the edge
         } 
      } 
   }
   
   # return the updated A, alpha, beta, delta, gamma, logit(pi), and log(mu)
   return(list(A       = A,
               alpha   = alpha,
               beta    = beta,
               delta   = delta,
               gamma   = gamma,
               logitPi = logitPi,
               logMu   = logMu))
}

# update tau via Gibbs sampling
update_tau = function(alpha, beta, delta, gamma, A, priors, temp = 1)
{
   p   = nrow(A)
   tau = rep(NA, 4)
   
   # sample tau_alpha
   tau[1] = rgamma(1, shape = (priors$tau[1] + 0.5 * sum(A) + temp - 1) / temp, rate = (priors$tau[2] + 0.5 * sum(alpha * alpha)) / temp)
   
   # sample tau_beta
   tau[2] = rgamma(1, shape = (priors$tau[1] + 0.5 * sum(A) + temp - 1) / temp, rate = (priors$tau[2] + 0.5 * sum(beta * beta)) / temp)
   
   # sample tau_delta
   tau[3] = rgamma(1, shape = (priors$tau[1] + 0.5 * p + temp - 1) / temp, rate = (priors$tau[2] + 0.5 * sum(delta * delta)) / temp)
   
   # sample tau_gamma
   tau[4] = rgamma(1, shape = (priors$tau[1] + 0.5 * p + temp - 1) / temp, rate = (priors$tau[2] + 0.5 * sum(gamma * gamma)) / temp)
   
   # return the sampled tau
   return(tau)
}

# update rho via Gibbs sampling
update_rho = function(A, priors, temp = 1)
{
   p = nrow(A)
   
   # sample rho
   rho = rbeta(1, shape1 = (priors$rho[1] + sum(A) + temp - 1) / temp, shape2 = (priors$rho[2] + p * (p - 1) - sum(A) + temp - 1) / temp)
   
   # return the sampled rho
   return(rho)
}

# evaluate log-likelihood of each observation for the j-th component of ZINBBN model
llik_ZINBBN_j = function(x, logitPi, logMu, psi)
{
   # calculate pi and lambda
   pi = exp(logitPi) / (1 + exp(logitPi))
   mu = exp(logMu)
   pi[is.nan(pi)] = 1
   mu[mu == Inf] = .Machine$double.xmax
   
   # evaluate and return log-likelihood of each observation
   llik = rep(0, length(x))
   llik[x == 0] = log(pi[x == 0] + (1 - pi[x == 0]) * dnbinom(0, size = 1 / psi, mu = mu[x == 0], log = FALSE))
   llik[x > 0]  = log(1 - pi[x > 0]) + dnbinom(x[x > 0], size = 1 / psi, mu = mu[x > 0], log = TRUE)
   return(llik)
}
