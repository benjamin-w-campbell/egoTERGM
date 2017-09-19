#' Estimation of ego-Temporal Exponential Random Graph Model (ego-TERGM) using Expectation Maximization (EM).
#'
#' This function estimates an ego-TERGM.
#' @param net The longitudinally observed network that an ego-TERGM will be fit on.  Must be presented as a list of networks.  Any vertex attributes should be attached to networks.
#' @param core_size The order of alters to include. The defaul value of one implies only looking at an ego's alters and the connections among them.
#' @param min_size  The minimum number of nodes an ego-network must achieve to be included.  Defaults to five.
#' @param roles The number of roles that should be fit.  Defaults to 3.
#' @param form The formula comprised of ERGM or TERGM terms used to distuinguish between clusters assignments.  Specified as a vector of comma seperated terms. No default.
#' @param add_drop Do nodes drop out of the network or enter it? If so, specify as the default TRUE.
#' @param directed Should the longitudinal network be treated as directed? If so, specify as the default TRUE.
#' @param edge_covariates Are edge covariates included in the form term? IF so, specify as TRUE.  No default.
#' @param seed The seed set to replicate analysis for pseudorandom number generator.
#' @param R The number of bootstrap replications that should be used for the estimation of a bootstrapped MPLE estimated TERGM for model initialization.  Defaults to 10.
#' @param parallel How the BTERGM should be computed, either with parallel processing ("multicore", "snow") or no parallel processing ("no").  Defaults to no.
#' @param ncpus The number of CPUs that should should be used for estimation, defaults to 1.
#' @param steps The number of default EM steps that should be taken, defaults to 50.
#' @param tol The difference in parameter estimates between EM iterations to determine if the algorithm has converged.  Defaults to 1e-6.
#' @return A list of model results, including lambda (the probability of assignments), group.theta (the roles by terms cluster centroids),
#'         EE.BIC (the Salter-Townshend and Murphy BIC cross-sectional BIC), TS.BIC (the Campbell BIC penalizing for time-steps),
#'        role_assignments (a data frame of the most likely assignments), reduced_networks (A list of the networks with excluded egos), and
#'        ego_nets (a list of ego-networks).
#' @keywords ego-TERGM
#' @references Campbell, Benjamin W. 2017. Inferring Latent Roles in Longitudinal Networks using the Ego-TERGM. Working Paper.
#' @examples
#' \dontrun{
#' # Code from xergm.common and their preparation of the Knecht network
#' library(xergm.common)
#' set.seed(1)
#'
#' data("knecht")
#'
#' for (i in 1:length(friendship)) {
#'  rownames(friendship[[i]]) <- paste("Student.", 1:nrow(friendship[[i]]), sep="")
#'  colnames(friendship[[i]]) <- paste("Student.", 1:nrow(friendship[[i]]), sep="")
#' }
#' rownames(primary) <- rownames(friendship[[1]])
#' colnames(primary) <- colnames(friendship[[1]])
#' sex <- demographics$sex
#' names(sex) <- rownames(friendship[[1]])
#' # step 2: imputation of NAs and removal of absent nodes:
#' friendship <- xergm.common::handleMissings(friendship, na = 10, method = "remove")
#' friendship <- xergm.common::handleMissings(friendship, na = NA, method = "fillmode")
#' # step 3: add nodal covariates to the networks
#' for (i in 1:length(friendship)) {
#'   s <- xergm.common::adjust(sex, friendship[[i]])
#'   friendship[[i]] <- network::network(friendship[[i]])
#'   friendship[[i]] <- network::set.vertex.attribute(friendship[[i]], "sex", s)
#'   idegsqrt <- sqrt(sna::degree(friendship[[i]], cmode = "indegree"))
#'   friendship[[i]] <- network::set.vertex.attribute(friendship[[i]],
#'                                                    "idegsqrt", idegsqrt)
#'   odegsqrt <- sqrt(sna::degree(friendship[[i]], cmode = "outdegree"))
#'   friendship[[i]] <- network::set.vertex.attribute(friendship[[i]],
#'                                                    "odegsqrt", odegsqrt)
#' }
#' sapply(friendship, network::network.size)
#' net <- friendship
#' rm(list=setdiff(ls(), "net"))
#'
#' ego_tergm_fit <- ego_tergm(net = net,
#'                           form = c("edges", "mutual", "triangle",
#'                                    "nodeicov('idegsqrt')", "nodeocov('odegsqrt')",
#'                                    "nodematch('sex')"),
#'                           core_size = 1,
#'                           min_size = 5,
#'                           roles = 3,
#'                           add_drop = TRUE,
#'                           directed = TRUE,
#'                           edge_covariates = FALSE,
#'                           seed = 12345,
#'                           R = 10,
#'                           parallel = "no",
#'                           ncpus = 1,
#'                           steps = 50,
#'                           tol = 1e-06)
#' }
#' @export

ego_tergm <- function(net = NULL,
                      form = NULL,
                      core_size = 1, min_size = 5, roles = 3, add_drop = TRUE, directed = TRUE, edge_covariates = FALSE,
                      seed = 12345,
                      R = 10, parallel = "no", ncpus = 1,
                      steps = 50, tol = 1e-6){

  # tested prior to 9/3/17

  N = max(unlist(lapply(net, function(x) network::network.size(x))))

  ########################################################################
  ### Start Functions
  ###  (for network prep)
  ########################################################################
  time_steps <- length(net)

  # extract the all vertices never appearing in a particular network but some network and add them as isolates
  if(add_drop==TRUE){
    vertices <- unique(unlist(lapply(net, function(x) network::get.vertex.attribute(x, 'vertex.names'))))

    # x = net
    add_setdiff <- function(x){
      ids   <- network::get.vertex.attribute(x, 'vertex.names')
      to_add <- setdiff(vertices,ids)
      new_vertex_ids <- c(ids, to_add)

      x <- network::add.vertices(x, length(to_add))
      x <- network::set.vertex.attribute(x, 'vertex.names', new_vertex_ids)

    }

    net <- lapply(net, function(x) add_setdiff(x))


    N = max(unlist(lapply(net, function(x) network::network.size(x))))
  } else {
    vertices <- unique(unlist(lapply(net, function(x) network::get.vertex.attribute(x, 'vertex.names'))))
  }


  YT <- lapply(net, function(x) network::as.matrix.network(x))

  ### find each ego-network; use K steps out from each node
  # use Y+t(Y) to jump backwards across links at the second step out

  if(directed == FALSE){
    YT2 <- lapply(YT, function(x) x+t(x))
  }

  if(directed == TRUE){
    YT2 <- YT

  }


  net2 <- lapply(YT2, function(x) network::network(x, directed = directed)) #network object based upon the network matrix y which takes y and transforms it by	causing nodes to "jump backwards across links at the second step"

  neighborhood_extract <- function(x){
    net_x <- sna::gapply(x,c(1,2),1:N,"*",1,distance=core_size)
  }

  xt<-lapply(net2, function(x) neighborhood_extract(x))
  #This abov line returns a vector obtained by applying a function to vertex neighborhoods ina  certain order. Here, this has the net2 object as the input graph, the margin of x to be used in calculating the network so that it is includin both rows and columns (total) = that's the second object in this function.  This does it for all nodes 1 through N, applying the star function with a maximum distance of K in which neighborhoods are taken.
  #The following for loop says for every node. apply a function over the list of vectors. I think this is just a way to change the matrix into a more readable x,y form
  for(ti in 1:length(net)){
    time_period <- xt[[ti]]
    x <- list()
    Y <- YT2[[ti]]
    for (n in 1:N){
      x[[n]]<-as.matrix(Y[c(n,time_period[[n]]),c(n,time_period[[n]])])
    }
    x<-lapply(x, function(x) network::as.network.matrix(x, directed=directed))
    xt[[ti]] <- x
  }


  # get this in xt[[ego]][[time-observation]] format
  ego_list <- list()
  for(n in 1:N){
    time_list <- list()
    for(ti in 1:length(net)){
      nit <- xt[[ti]][[n]]
      time_list[[ti]] <- nit
    }
    ego_list[[n]] <- time_list
  }
  xt <- ego_list

  # only look at egos bigger than MINSIZE.  Here I think the MINSIZE is however many connections they have
  if (min_size>1)
  {
    keep_mat <- matrix(ncol = length(net), nrow = N)
    for(i in 1:length(xt)){
      ego <- xt[[i]]
      keep <- lapply(ego, network::network.size)>=min_size
      keep_mat[i,] <- keep
    }

    # all_net<-net
    red_net_list <- list()
    if(length(names(net)) > 0){
      net <- unname(net)
    }

    for(t in 1:length(net)){
      time_slice <- net[[t]]
      red_net <- network::as.network(as.matrix(network::as.sociomatrix(time_slice)[(1:N)[keep_mat[,t]],(1:N)[keep_mat[,t]]]), directed=directed)
      for(att in network::list.vertex.attributes(time_slice)){
        network::set.vertex.attribute(red_net,att,network::get.vertex.attribute(time_slice,att)[keep_mat[,t]])
      }
      if(edge_covariates == TRUE){
        for(att_e in network::list.edge.attributes(time_slice)){
          # setting edge attributes harder given how they're stored
          el <- as.matrix(time_slice,matrix.type="edgelist")
          el[,1] <- network::get.vertex.attribute(time_slice, 'vertex.names')[el[,1]]
          el[,2] <- network::get.vertex.attribute(time_slice, 'vertex.names')[el[,2]]
          el <- cbind(el, network::get.edge.attribute(time_slice, att_e))

          el_red <- as.matrix(red_net, matrix.type="edgelist")
          el_red[,1] <- network::get.vertex.attribute(red_net, 'vertex.names')[el_red[,1]]
          el_red[,2] <- network::get.vertex.attribute(red_net, 'vertex.names')[el_red[,2]]
          el_red <- cbind(el_red, NA)

          # check to make sure that el's are same order
          if(all(el_red[,1] %in% el[,2])){
            el_red <- cbind(el_red[,2], el_red[,1], NA)
          }

          # aad a third column to el_red which contains the edge attribute from el when edges were kept
          m3 <- rbind(el, el_red)
          el <- m3[duplicated(m3[,1:2], fromLast = TRUE), , drop = TRUE]

          network::set.edge.attribute(red_net,att_e,el[,3])

        }
      }

      red_net_list[[t]] <- red_net
    }
    reduced_networks <- red_net_list

    net <- red_net_list

    # reduce xt list by keep mat
    for(i in 1:nrow(keep_mat)){
      ego <- xt[[i]]
      ego <- ego[keep_mat[i,]]
      xt[[i]] <- ego
    }
    N=length(xt)

    for(i in 1:length(xt)){
      nets <- xt[[i]]
      for(t in 1:length(nets)){
        if(length(nets) > 1){
          vertex_ids <- network::get.vertex.attribute(xt[[i]][[t]], 'vertex.names')
          indices <- which(vertex_ids %in% network::get.vertex.attribute(net[[t]], 'vertex.names'))
          for(att in network::list.vertex.attributes(net[[t]])){
            network::set.vertex.attribute(xt[[i]][[t]], att, network::get.vertex.attribute(net[[t]], att)[indices])
          }
          if(edge_covariates == TRUE){
            for(att_e in network::list.edge.attributes(net[[t]])){
              el <- as.matrix(net[[t]],matrix.type="edgelist")
              el[,1] <- network::get.vertex.attribute(net[[t]], 'vertex.names')[el[,1]]
              el[,2] <- network::get.vertex.attribute(net[[t]], 'vertex.names')[el[,2]]
              el <- cbind(el, network::get.edge.attribute(net[[t]], att_e))

              el_red <- as.matrix(xt[[i]][[t]], matrix.type="edgelist")
              el_red[,1] <- network::get.vertex.attribute(xt[[i]][[t]], 'vertex.names')[el_red[,1]]
              el_red[,2] <- network::get.vertex.attribute(xt[[i]][[t]], 'vertex.names')[el_red[,2]]
              el_red <- cbind(el_red, NA)

              if(all(el_red[,1] %in% el[,2])){
                el_red <- cbind(el_red[,2], el_red[,1], NA)
              }

              if(network::network.size(xt[[i]][[t]]) == 1){
                el_red <- matrix(nrow = 0, ncol = 3)
              }

              m3 <- rbind(el, el_red)
              el <- m3[duplicated(m3[,1:2], fromLast = TRUE), , drop = TRUE]

              if(class(el) == "numeric" & length(el) == 3){
                el <- t(as.matrix(el))
              }

              network::set.edge.attribute(xt[[i]][[t]],att_e,el[,3])
            }

          }

        }
        if(length(nets) == 1 & t == 1){
          vertex_ids <- network::get.vertex.attribute(nets[[1]], 'vertex.names')
          indices <- which(network::get.vertex.attribute(net[[t]], 'vertex.names') %in% vertex_ids) #
          for(att in network::list.vertex.attributes(net[[t]])){
            network::set.vertex.attribute(nets[[1]], att, network::get.vertex.attribute(net[[t]], att)[indices])
          }
          if(edge_covariates == TRUE){
            for(att_e in network::list.edge.attributes(net[[t]])){
              el <- as.matrix(net[[t]],matrix.type="edgelist")
              el[,1] <- network::get.vertex.attribute(net[[t]], 'vertex.names')[el[,1]]
              el[,2] <- network::get.vertex.attribute(net[[t]], 'vertex.names')[el[,2]]
              el <- cbind(el, network::get.edge.attribute(net[[t]], att_e))

              el_red <- as.matrix(nets[[1]], matrix.type="edgelist")
              el_red[,1] <- network::get.vertex.attribute(nets[[1]], 'vertex.names')[el_red[,1]]
              el_red[,2] <- network::get.vertex.attribute(nets[[1]], 'vertex.names')[el_red[,2]]
              el_red <- cbind(el_red, NA)

              if(all(el_red[,1] %in% el[,2])){
                el_red <- cbind(el_red[,2], el_red[,1], NA)
              }

              if(network::network.size(xt[[i]][[t]]) == 1){
                el_red <- matrix(nrow = 0, ncol = 3)
              }

              m3 <- rbind(el, el_red)
              el <- m3[duplicated(m3[,1:2], fromLast = TRUE), , drop = TRUE]

              if(class(el) == "numeric" & length(el) == 3){
                el <- t(as.matrix(el))
              }

              network::set.edge.attribute(nets[[1]],att_e,el[,3])
            }

          }

          xt[[i]] <- nets
        }
        if(length(nets) == 0){
          xt[[i]] <- NA
        }
      }
    }
    null_networks <- which(is.na(xt))
    dropped <- vertices[null_networks]
    print(paste(dropped, "removed from x list due to no identifiable networks"))
    remaining_vertices <- vertices[-c(null_networks)]
    xt <- xt[!is.na(xt)]

    xt <- xt[!is.na(xt)]


    rm(keep_mat, ego, keep, time_slice, red_net, i, t, red_net_list, nets, vertex_ids, indices)
  } else {
    stop("Minimum size must be larger than 1.  Please adjust the min_size argument accordingly.")
  }
  x <- xt
  rm(xt)

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
  EE.BIC<-out$EE.BIC
  TS.BIC <- out$TS.EE.BIC
  z<-apply(lambda, 1, which.max)
  roles_out <- data.frame(Id = remaining_vertices,
                          Role = z)

  return(list(model.fit = "egoTERGM",
              lambda = lambda,
              group.theta = group.theta,
              EE.BIC = EE.BIC,
              TS.BIC = TS.BIC,
              role_assignments = roles_out,
              reduced_networks = reduced_networks,
              form = form,
              ego_nets = x))

}
