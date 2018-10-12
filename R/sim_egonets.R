#' Simulation of longitudinal ego-networks according to a mixture of data generating processes and fitting of ego-TERGM to that network.
#'
#' This function simulated longitudinal ego-networks and estimates an ego-TERGM. Useful for monte carlos and proofs of concept.
#' @param form A vector of ERGM terms that should be used to generate the networks.
#' @param params A "roles" by "form" matrix of ERGM simulation parameters with groups defined on the row and model terms defined on the column.
#' @param roles An integer for the number of distinct mixture groups that should be simulated.
#' @param N_per_role An integer for the number of different longitudinally observed ego-networks that should be simulated per role in roles.
#' @param t_steps An integer for the number of time steps that each ego-network should be observed across.
#' @param egonet_size An integer for the size of each ego-network simulated.
#' @param seed The seed set to replicate analysis for pseudorandom number generator.
#' @param R The number of bootstrap replications that should be used for the estimation of a bootstrapped MPLE estimated TERGM for model initialization.  Defaults to 10.
#' @param forking If parallelization via forking should be used (TRUE) or if no parallel processing should be used (FALSE).  Currently, sockets are not supported.
#' @param ncpus The number of CPUs that should should be used for estimation, defaults to 1.
#' @param steps The number of default EM steps that should be taken, defaults to 50.
#' @param tol The difference in parameter estimates between EM iterations to determine if the algorithm has converged.  Defaults to 1e-6.
#' @return A list of simulated ego-networks and the output of the ego_tergm function fit to this.
#' @keywords simulation networks
#' @references{
#' #' Campbell, Benjamin W. (2018):
#'  Inferring Latent Roles in Longitudinal Networks.
#'  \emph{Political Analysis} 26(3): 292-311.  \url{https://doi.org/10.1017/pan.2018.20}
#'
#'  Salter-Townshend, Michael and Thomas Brendan Murphy. (2015):
#'  Role Analysis in Networks using Mixtures of Exponential Random Graph Models.
#'  \emph{Journal of Computational and Graphical Statistics} 24(2): 520-538. \url{https://doi.org/10.1080/10618600.2014.923777}
#' }
#' @examples
#' \donttest{
#' net <- sim_egonets(form = c("edges", "gwesp(0.8,fixed=TRUE)", "gwdegree(decay=0.8,fixed=TRUE)"),
#'                   params = rbind(c(-3,1,0), c(-1,-2,-1), c(-2,0,2)),
#'                   roles = 3,
#'                   N_per_role = 10,
#'                   t_steps = 3,
#'                   egonet_size = 20,
#'                   seed = 12345,
#'                   R = 30,
#'                   forking = FALSE,
#'                   ncpus = 1,
#'                   steps = 50,
#'                   tol = 1e-6)
#' }
#' @export

sim_egonets <- function(form = NULL, params = NULL, roles = NULL, N_per_role = NULL, t_steps = NULL, egonet_size = NULL,
                        R = 10, forking = FALSE, ncpus = 1, steps = 50, tol = 1e-6, seed = 12345){
  cat("Start Time:", format(Sys.time(), "%a %b %d %X %Y"), "\n")
  set.seed(seed)
  if(length(form) != ncol(params)){
    stop("The form argument is not of same length as the number of columns in the params matrix.  Please provide one column of simulation parameters for each term in form")
  }

  if(roles != nrow(params)){
    stop("The number of roles provided is not equal to the number of rows in the params matrix.  Please provide one row of simulation parameters for each role")
  }

  N <- roles*N_per_role

  time_steps <- t_steps

  cat("Simulating ego-networks according to parameters provided.", "\n")
  ### assign the egos to groups
  sim.K<-rep(NaN,N)
  sim.list <- list()
  for (r in 1:roles){
    sim.K[((r-1)*N_per_role+1):(r*N_per_role)] <- rep(r, N_per_role)
  }

  for(i in 1:length(sim.K)){
    r <- sim.K[i]
    sim.list[[i]] <- rep(r, t_steps)
  }

  # sim ego nets
  time_list <- list()
  sim_egos <- function(k){
    sim.K <- sim.list[[k]]
    for(i in 1:length(sim.K)){
      net <- network::network.initialize(n = egonet_size, directed = FALSE)
      time_list[[i]]<- simulate(object = as.formula(paste("net", "~", paste(form, collapse="+"), sep="")), coef = c(params[(1:roles)[sim.K[i]],]))
    }
    return(time_list)
  }

  if(forking == TRUE){
    sim.K <- parallel::mclapply(seq_along(sim.K), sim_egos, mc.cores = ncpus)
  } else {
    sim.K <- lapply(seq_along(sim.K), sim_egos)
  }

  ########################################################################
  ### Likelihood functions
  ########################################################################

  pseudo.loglikelihood<-function(S,tmp.theta)  # Establish a pseudo-loglikelihood function
  {
    loglike=sapply(S, function (S) sum(stats::dbinom(S$response*S$weights,S$weights,
                                                     boot::inv.logit(S$offset+as.matrix(S$predictor)%*%tmp.theta),log=TRUE),na.rm=1))
    loglike=sum(loglike)
    if(!is.finite(loglike)|loglike<LOWESTLL)
      loglike=LOWESTLL# avoids numerical errors
    return(loglike)
  }
  # returns the negative so that optimization is towards maximum
  n.pseudo.loglikelihood<-function(S,tmp.theta){
    -pseudo.loglikelihood(S,tmp.theta)
  }

  approx.loglikelihood<-function(S,tmp.theta,old.theta,M,form,ll0){
    pseudo.loglikelihood(S,tmp.theta)
  }
  # loglikelihood summed across all groups
  Mstepfunction<-function(tmp.theta,S,N,lambda,TAU,roles)
  {
    tmp.theta<-matrix(tmp.theta,nrow=roles)
    ans=0
    for (g in 1:nrow(tmp.theta))
      ans = ans + mstepfunction(tmp.theta[g,],S,N,lambda[,g],TAU[g])
    return(ans)
  }
  n.approx.loglikelihood<-function(S,tmp.theta,old.theta,M,form,ll0){
    -approx.loglikelihood(S,tmp.theta,old.theta,M,form,ll0)
  }


  mstepfunction<-function(tmp.theta,S,N,lambda,TAU,old.theta,M,form,ll0){
    sum(lambda * (log(TAU) + sapply(S,approx.loglikelihood,tmp.theta,old.theta,M,form,ll0)))
  }

  n.mstepfunction<-function(tmp.theta,S,N,lambda,TAU,old.theta,M,form,ll0){
    -mstepfunction(tmp.theta,S,N,lambda,TAU,old.theta,M,form,ll0)
  }

  x = sim.K



  ########################################################################
  ### Initialization functions
  ########################################################################
  init.egoergm <- function(form = NULL, roles = roles, p.ego.terms=NULL, R=R, seed = seed, N = length(x)){
    set.seed(seed)
    Nterms<-length(form) #specifying a "nterms" object as the length of the number of ego.terms

    #This specifies object "Form" as an object that is a new formula that is a function of an indexed xi, according to
    #the ergm formula that is specied above

    ######### fit all ego-networks #############
    #if p.ego terms is null or empty, then do the following...
    if (is.null(p.ego.terms)){
      # new
      form = form
      fit_btergm_local <- function(i, form = NULL){
        form <- stats::formula(paste("x[[i]]~", paste(form,collapse="+"),sep=""))
        #Specify object ergmformula - paste command has three arguments - the stuff you want to paste together
        #, sep is how to separate them, and collapse is if you want smush it together.  This specifies ergmformula
        #as ~ pasted together with ego terms, separated by " "
        # form<-stats::update.formula(paste("x[[i]]",ergmformula),x[[i]] ~ .)
        fit = btergm::btergm(formula = form, R=R, verbose = FALSE)@coef
        return(list(fit))
      }
      if(forking == TRUE){
        theta <- parallel::mclapply(seq_along(x), function(i) fit_btergm_local(i, form = form), mc.cores = ncpus)
      } else {
        theta <- lapply(seq_along(x), function(i) fit_btergm_local(i, form = form))
      }

      theta<- do.call(rbind, unlist(theta, recursive = FALSE))
      # next couple of lines very ad-hoc but not an issue post EM convergence.
      theta[is.na(theta)]<-0
      theta[theta==-Inf]<- -1e6
      theta[theta==Inf]<- 1e6
      ############# initialisation ##############
      # k-means Clustering start #if statement - if G>1, then the following is true and initial.clusters is set equal to
      #the kmeans of which takes the theta matrixm, for the number of centers, it is a function of theta, 2, and some other stuff?
      #I don't really understand this portion too much - go back to the article and see where this comes in?
      if(roles>1){
        # initial cluster assignments look good
        initial.clusters<-stats::kmeans(theta, roles, centers=apply(theta, 2, tapply,stats::cutree(stats::hclust(stats::dist(theta)),roles), mean),nstart=5)
        group.theta<-initial.clusters$centers
        # initial centers don't look great
        group.theta<-matrix(group.theta,roles,Nterms)
      }
      if(roles==1){
        group.theta<-matrix(apply(theta,2,mean),nrow=1)
      }
      lambda<-matrix(NaN,N,roles)
      for (i in 1:N){
        for (g in 1:roles){
          lambda[i,g]<-1/(sqrt(t(group.theta[g,]-theta[i,])%*%(group.theta[g,]-theta[i,]))+tol) # nugget to avoid zero
        }
      }
      lambda<-lambda/apply(lambda,1,sum) # normalise lambda
      cat("Finished kmeans initialization.", "\n")
    }
    return(list(theta=theta, group.theta=group.theta, lambda=lambda, roles=roles))
  }


  init<-init.egoergm(form = form, roles = roles, p.ego.terms = NULL, R = R, seed = seed)

  ########################################################################
  ### Calculating and save change statistics
  ########################################################################

  cat("Calculating all network statistics.", "\n")

  calculate_change_stats <- function(i){
    ego_lists <- list()
    temp <- x[[i]]
    for(i in 1:length(temp)) {
      net<-temp[[i]]
      cs <- ergm::ergmMPLE(stats::as.formula(paste("net",ergmformula)))
      cs$offset <- -log(network::network.size(net))
      ego_lists[[i]] <- cs
    }
    return(ego_lists)
  }

  ergmformula <- paste("~", paste(form,collapse="+"),sep="") # Establish function ergm formula that includes the ego.terms object

  if(forking == TRUE){
    obs.S <- parallel::mclapply(seq_along(x), calculate_change_stats, mc.cores = ncpus)
  } else {
    obs.S <- lapply(seq_along(x), calculate_change_stats)
  }

  cat("Network statistics calculated.", "\n")

  ########################################################################
  ### EM Algorithm
  ########################################################################

  fit.mix.egoergm<-function(form, init, obs.S, roles, p.ego.terms=NULL, p.group.theta=NULL, N = length(obs.S), x = x, STEPS = steps, M)
    # Specify a function to perform EM algorithm to find ego-ERGM parameters and memberships.
  {
    Nterms<-length(form)
    ergmformula <- paste("~", paste(form,collapse="+"),sep="")
    form<-stats::update.formula(stats::as.formula(paste("x[[i]]",ergmformula)),x[[i]] ~ .)
    lambda<-init$lambda
    group.theta<-init$group.theta
    TAU<-apply(lambda,2,sum)/N
    LL<-NaN
    cat("iterations\tlog-likelood\n")
    ############# E-M iterations ##############
    for (l in 1:STEPS)
    {
      oldLL<-LL
      #cat(l, " ", sep="")
      # E-step
      #cat("E-step", l, "\t")
      old_lambda<-lambda
      for (i in 1:N){
        for (g in 1:roles){
          # so, g<-2 then the ll is artificially high
          lambda[i,g]<- log(TAU[g]) + approx.loglikelihood(obs.S[[i]],group.theta[g,],group.theta[g,]*0,M,form,0)
        }
        # need to fix approx.log
        lambda[i,]<-lambda[i,]-max(lambda[i,])
      }
      lambda<-exp(lambda)
      lambda<-lambda/apply(lambda,1,sum)
      lambda[is.na(lambda)]<-1/roles
      tmp<-apply(lambda,1,sum)
      lambda<-lambda/tmp # normalise lambda
      lambda[lambda==0]<-1e-8; lambda<-lambda/tmp
      TAU<-apply(lambda,2,sum)/N
      # M-step
      #cat("M-step", l, "\t")
      LL<-Mstepfunction(group.theta, obs.S, N, lambda, TAU, roles)
      cat(l, "\t", sprintf("%.10f",LL),"\n")
      old_group.theta<-group.theta
      for (g in 1:roles)
        group.theta[g,]<-stats::nlm(n.mstepfunction,group.theta[g,],S=obs.S,N=N,lambda=lambda[,g],TAU=TAU[g],group.theta[g,]*0,M,form,0)$estimate
      #if (max(abs(group.theta-old_group.theta))<tol)
      if (l>1)
        if ((LL-oldLL) < tol)
        {
          cat("Converged after ", l, "iterations\n")
          break
        }
    }
    #LL<-Mstepfunction(group.theta,obs.S,N,lambda,TAU,G)
    # higher is better
    CS.EE.BIC<- 2*LL - (2*roles*Nterms+roles-1)*log(N) # usual BIC
    TS.EE.BIC<- 2*LL - (2*roles*Nterms+roles-1)*log(N*time_steps)
    return(list(EE.BIC=CS.EE.BIC, TS.EE.BIC=TS.EE.BIC, theta=group.theta, lambda=lambda, roles=roles))
  }

  LOWESTLL=-1e8

  cat("EM algorithm starting.", "\n")

  out <- fit.mix.egoergm(form=form,init=init,obs.S,roles) # run the EM algorithm
  lambda<-out$lambda
  group.theta<-out$theta
  EE.BIC<-out$EE.BIC
  TS.BIC <- out$TS.EE.BIC
  z<-apply(lambda, 1, which.max)
  roles_out <- data.frame(Id = paste0("v", 1:N),
                          Role = z)

  cat("EM algorithm completed.", "\n")
  cat("Done.", "\n")
  cat("Completed Time:", format(Sys.time(), "%a %b %d %X %Y"), "\n")

  return(list(lambda = lambda,
              group.theta = group.theta,
              EE.BIC = EE.BIC,
              TS.BIC = TS.BIC,
              role_assignments = roles_out,
              reduced_networks = sim.K,
              form = form))
}
