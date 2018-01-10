#' Simulation of overlapping longitudinal ego-networks according to a mixture of data generating processes and fitting of ego-TERGM to that network.
#'
#' This function simulated longitudinal ego-networks that are dependent and overlap, and estimates an ego-TERGM. This is done according to the conventional DGP used in the manuscript.  Useful for monte carlos an overlap. Note, to simplify code, the DGP and number of roles are fixed for this example. 5 networks overlap per directed DGP combination (so, 30 total)
#' @param N_per_role An integer for the number of different longitudinally observed ego-networks that should be simulated per role in roles.
#' @param t_steps An integer for the number of time steps that each ego-network should be observed across.
#' @param egonet_size An integer for the size of each ego-network simulated.
#' @param seed The seed set to replicate analysis for pseudorandom number generator.
#' @param prop_fixed A proportion of ties in overlapping networks that should be shared between data generating processes.
#' @param R The number of bootstrap replications that should be used for the estimation of a bootstrapped MPLE estimated TERGM for model initialization.  Defaults to 10.
#' @param parallel How the BTERGM should be computed, either with parallel processing ("multicore", "snow") or no parallel processing ("no").  Defaults to no.
#' @param ncpus The number of CPUs that should should be used for estimation, defaults to 1.
#' @param steps The number of default EM steps that should be taken, defaults to 50.
#' @param tol The difference in parameter estimates between EM iterations to determine if the algorithm has converged.  Defaults to 1e-6.
#' @param seed The seed for simulation
#' @return A list of simulated ego-networks and the output of the ego_tergm function fit to this.
#' @keywords simulation networks
#' @keywords internal
#'# @references Campbell, Benjamin W. 2017. Inferring Latent Roles in Longitudinal Networks using the Ego-TERGM. Working Paper.
#' @examples
#' \dontrun{
#' net <- sim_overlapping_egonets(N_per_role = 10,
#'                   t_steps = 3,
#'                   egonet_size = 20,
#'                   R = 30,
#'                   parallel = "multicore",
#'                   ncpus = 3,
#'                   steps = 50,
#'                   tol = 1e-6)
#' }
sim_overlap_egonets <- function(N_per_role = NULL, t_steps = NULL, egonet_size = NULL,
                        R = 10, parallel = "no", ncpus = 1, steps = 50, tol = 1e-6, seed = 12345){

  form = c("edges", "gwesp(0.8,fixed=TRUE)", "gwdegree(decay=0.8,fixed=TRUE)")
  params = rbind(c(-3,1,0), c(-1,-2,-1), c(-2,0,2))
  roles = 3
  set.seed(seed)

  if(prop_fixed > 1 | prop_fixed < 0){
    stop("The prop_fixed argument must be a proportion between 0 and 1")
  }

  if(length(form) != ncol(params)){
    stop("The form argument is not of same length as the number of columns in the params matrix.  Please provide one column of simulation parameters for each term in form")
  }

  if(roles != nrow(params)){
    stop("The number of roles provided is not equal to the number of rows in the params matrix.  Please provide one row of simulation parameters for each role")
  }

  N <- roles*N_per_role

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

  ### simulate the ego
  NN = egonet_size # number of nodes per ego
  ergmformula <- paste("~", paste(form, collapse="+"), sep="")
  sim.net<-list()
  sim.x <- list()
  for(k in 1:length(sim.K)){
    sim.K <- sim.list[[k]]
    for(i in 1:t_steps){
      temp <- simulate.formula(as.formula(paste("network(NN,directed=FALSE)", ergmformula)), coef = c(params[(1:roles)[sim.K[i]],]))
      # add node labels (will help with creating overlap)
      temp <- set.vertex.attribute(temp, 'vertex.names', paste0("node", 1:30))
      sim.net[[i]] <- temp
    }
    sim.x[[k]] <- sim.net
  }
  sim.K <- sim.x

  # Now program in overlap, this is done by taking N_overlap number of networks from each combination of mixing groups
    # and then modifying some portion of their edges to be fixed.  Then homophily values are assigned upon members to those fixed ties.

  combinations <- expand.grid(1:roles, 1:roles)
  combinations <- combinations[combinations$Var1 < combinations$Var2,]

  # empty edgelist
  el <- expand.grid(paste0("node", 1:30), paste0("node", 1:30))
  homophily_df <- data.frame(id = paste0("node", 1:30))
  # now remove loops
  el <- el[el$Var1 != el$Var2,]
  # make directed
  el <- el[as.numeric(gsub("node", "", el$Var1))<as.numeric(gsub("node", "", el$Var2)),]

  # function to do everything
  make_overlap <- function(net, fixed_ties, homophily_df){
    net_temp <- net
    net_temp[as.matrix(fixed_ties[,c(1:2)])] <- 1

    net <- set.vertex.attribute(net_temp, 'homophily', homophily_df$homophilous)

    return(net)
  }

  # First set of overlapping egonetworks
  indices_1to2 <- 26:35
  el_1 <- el
  el_1$fixed_tie <- rbinom(nrow(el_1), size = 1, prob = prop_fixed)
  el_1 <- el_1[el_1$fixed_tie == 1,]
  homophilous <- as.character(sort(unique(el_1$Var1, el_1$Var2)))
  homophily_df_1 <- homophily_df
  homophily_df_1$homophilous <- ifelse(homophily_df$id %in% homophilous, 1, 0)

  nets_1 <- sim.K[indices_1to2]
  for(i in 1:length(nets_1)){
    net_i <- nets_1[[i]]
    net_i <- lapply(net_i, function(x) make_overlap(net = x, fixed_ties = el_1, homophily_df = homophily_df_1))
    nets_1[[i]] <- net_i
  }
  sim.K[indices_1to2] <- nets_1

  # second set of overlapping egonetworks
  indices_2to3 <- 56:65
  el_2 <- el
  el_2$fixed_tie <- rbinom(nrow(el_2), size = 1, prob = prop_fixed)
  el_2 <- el_2[el_2$fixed_tie == 1,]
  homophilous <- as.character(sort(unique(el_2$Var1, el_2$Var2)))
  homophily_df_2 <- homophily_df
  homophily_df_2$homophilous <- ifelse(homophily_df$id %in% homophilous, 1, 0)

  nets_2 <- sim.K[indices_2to3]
  for(i in 1:length(nets_2)){
    net_i <- nets_2[[i]]
    net_i <- lapply(net_i, function(x) make_overlap(net = x, fixed_ties = el_2, homophily_df = homophily_df_2))
    nets_2[[i]] <- net_i
  }
  sim.K[indices_2to3] <- nets_2

  # third set of overlapping egonetworks
  indices_1to3 <- c(1:5, 86:90)
  el_3 <- el
  el_3$fixed_tie <- rbinom(nrow(el_3), size = 1, prob = prop_fixed)
  el_3 <- el_3[el_3$fixed_tie == 1,]
  homophilous <- as.character(sort(unique(el_3$Var1, el_3$Var2)))
  homophily_df_3 <- homophily_df
  homophily_df_3$homophilous <- ifelse(homophily_df$id %in% homophilous, 1, 0)

  nets_3 <- sim.K[indices_1to3]
  for(i in 1:length(nets_3)){
    net_i <- nets_3[[i]]
    net_i <- lapply(net_i, function(x) make_overlap(net = x, fixed_ties = el_3, homophily_df = homophily_df_3))
    nets_3[[i]] <- net_i
  }
  sim.K[indices_1to3] <- nets_3

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
  #Specify function in terms of ego.terms and G
  init.egoergm <- function(form = NULL, roles = roles, p.ego.terms=NULL, R=R, parallel = parallel, ncpus = ncpus, seed = seed, N = length(x)){
    set.seed(seed)
    Nterms<-length(form) #specifying a "nterms" object as the length of the number of ego.terms
    ergmformula <- paste("~", paste(form,collapse="+"),sep="")
    #Specify object ergmformula - paste command has three arguments - the stuff you want to paste together
    #, sep is how to separate them, and collapse is if you want smush it together.  This specifies ergmformula
    #as ~ pasted together with ego terms, separated by " "

    form<-ergm::ergm.update.formula(stats::as.formula(paste("x[[i]]",ergmformula)),x[[i]] ~ .)
    #This specifies object "Form" as an object that is a new formula that is a function of an indexed xi, according to
    #the ergm formula that is specied above

    ######### fit all ego-networks #############
    #if p.ego terms is null or empty, then do the following...
    if (is.null(p.ego.terms)){
      theta<-matrix(0,N,Nterms) #Theta is set equal to a matrix filled with zeros, with N rows and Nterms columns
      fit.list <- list()
      for (i in 1:length(x)){#For loop indexed from 1 to N - for every indexed i in the theta matrix, do ergmMPLE in the form object,
        #with the output of the fit, which is a coefficient for the model terms for each N.
        # for cross sectional ERGM, ergmMPLE would be sufficient, here we need BTERGM
        #theta[i,]<-ergmMPLE(form,output="fit")$coef

        fit <- btergm::btergm(form, R=R, parallel=parallel, ncpus=ncpus, verbose = FALSE)
        theta[i,] <- btergm::coef(fit)
        fit.list[[i]] <- fit
      }
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
      print("Finished kmeans initialisation")
    }
    return(list(theta=theta, group.theta=group.theta, lambda=lambda, roles=roles, fit.list=fit.list))
  }


  init<-init.egoergm(form = form, roles = roles, p.ego.terms = NULL, R = R, parallel = parallel, ncpus = ncpus, seed = seed)

  ########################################################################
  ### Calculating and save change statistics
  ########################################################################

  ergmformula <- paste("~", paste(form,collapse="+"),sep="") # Establish function ergm formula that includes the ego.terms object
  obs.S<-list() #obs.S list that will be useful for the below for loop
  print("Calculating all network statistics...")
  for (n in 1:length(x)) {# for loop started
    ego_lists <- list()
    temp <- x[[n]]
    for(i in 1:length(temp)) {
      net<-temp[[i]]
      cs <- ergm::ergmMPLE(stats::as.formula(paste("net",ergmformula)))
      cs$offset <- -log(network::network.size(net))
      ego_lists[[i]] <- cs
    }
    obs.S[[n]] <- ego_lists
  }

  time_steps <- t_steps

  ########################################################################
  ### EM Algorithm
  ########################################################################

  fit.mix.egoergm<-function(form, init, obs.S, roles, p.ego.terms=NULL, p.group.theta=NULL, N = length(obs.S), x = x, STEPS = steps, M)
    # Specify a function to perform EM algorithm to find ego-ERGM parameters and memberships.
  {
    Nterms<-length(form)
    ergmformula <- paste("~", paste(form,collapse="+"),sep="")
    form<-ergm.update.formula(stats::as.formula(paste("x[[i]]",ergmformula)),x[[i]] ~ .)
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

  out <- fit.mix.egoergm(form=form,init=init,obs.S,roles) # run the EM algorithm
  lambda<-out$lambda
  group.theta<-out$theta
  EE.BIC<-out$CS.EE.BIC
  TS.BIC <- out$TS.EE.BIC
  z<-apply(lambda, 1, which.max)
  roles_out <- data.frame(Id = paste0("v", 1:N),
                          Role = z)

  return(list(lambda = lambda,
              group.theta = group.theta,
              EE.BIC = EE.BIC,
              TS.BIC = TS.BIC,
              role_assignments = roles_out,
              reduced_networks = sim.K,
              form = form))

}
