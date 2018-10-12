#' Stability validation of ego-Temporal Exponential Random Graph Model (ego-TERGM) fit.
#'
#' This function examines the stability of ego-TERGM results to the pooling window.  It takes some proportion of the window's network and compares the result of a model fit on these time periods to the original fit.
#' @param ego_tergm_fit The output of a previously fit ego-TERGM fit using the ego_tergm function.  This is the model that stability validation will be performed on.
#' @param splitting_probability A value from 0 to 1 that determines the probability that any given network is assigned to the comparison group.
#' @param seed The seed set to replicate analysis for pseudorandom number generator.
#' @param R The number of bootstrap replications that should be used for the estimation of a bootstrapped MPLE estimated TERGM for model initialization.  Defaults to 10.
#' @param forking If parallelization via forking should be used (TRUE) or if no parallel processing should be used (FALSE).  Currently, sockets are not supported.
#' @param ncpus The number of CPUs that should should be used for estimation, defaults to 1.
#' @param steps The number of default EM steps that should be taken, defaults to 50.
#' @param tol The difference in parameter estimates between EM iterations to determine if the algorithm has converged.  Defaults to 1e-6.
#' @return Returns comparison_table (a matrix of cross-tabulation results to compare common cluster assignments or if incompatible a table of relative proportions sorted by value to allow for comparisons under set incompatibility and label switching), networks_sampled (which networks were included in the new validation sample),
#'         comparison_lambda (the matrix of role assignments for validation networks), comparison_group.theta (centroids for validation networks),
#'         comparison_EE.BIC (Salter-Townshend and Murphy BIC that doesn't penalize for longitudinal networks for validation networks),
#'         comparison_TS.BIC (BIC that penalizes for longitudinal networks for validation networks), comparison_role_assignments (role assignments for validation networks),
#'         and comparison_ego_nets (validation ego-networks).  Note that labels may switch.
#' @keywords ego-TERGM validation
#' @references{
#'  Campbell, Benjamin W. (2018):
#'  Inferring Latent Roles in Longitudinal Networks.
#'  \emph{Political Analysis} 26(3): 292-311.  \url{https://doi.org/10.1017/pan.2018.20}
#' }
#' @examples
#' \donttest{
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
#'                           forking = FALSE,
#'                           ncpus = 1,
#'                           steps = 50,
#'                           tol = 1e-06)
#'
#' stability_check <- stability_validation(ego_tergm_fit = ego_tergm_fit, seed = 614)
#' print(stability_check$comparison_table)
#' }
#' @export
#'
stability_validation <- function(ego_tergm_fit = NULL,
                                 splitting_probability = 0.5,
                                 seed = 12345,
                                 R = 10, forking = FALSE, ncpus = 1,
                                 steps = 50, tol = 1e-6){
  cat("Start Time:", format(Sys.time(), "%a %b %d %X %Y"), "\n")

  if(!is.null(seed)){
    set.seed(seed)
  }

  # create the random sample to fit a new ego-TERGM to:
  # First, create a set of indices to sample
  sample_indices <- rbinom(n = length(ego_tergm_fit$net), size = 1, prob = splitting_probability)
  if(sum(sample_indices) == 0 | sum(sample_indices) == length(sample_indices)){
    stop("Either no networks are resampled, or all networks are resampled.  Please either change the seed or or modify the splitting probability.")
  }

  # Second, create the subset
  subset_networks <- ego_tergm_fit$net[which(sample_indices == 1)]

  net <- subset_networks
  orig_nets <- subset_networks
  time_steps <- length(net)

  cat("Data formatting started.", "\n")
  # tested prior to 9/3/17

  N <-  max(unlist(lapply(net, function(x) network::network.size(x))))

  core_size = ego_tergm_fit$core_size
  min_size = ego_tergm_fit$min_size
  roles = ego_tergm_fit$roles
  add_drop = ego_tergm_fit$add_drop
  directed = ego_tergm_fit$directed
  edge_covariates = ego_tergm_fit$edge_covariates
  form = ego_tergm_fit$form

  if (min_size<=1){
    stop("Minimum size must be greater than 1.  Please adjust the min_size argument to a value greater than 1.")
  }
  ########################################################################
  ### Start Functions
  ###  (for network prep)
  ########################################################################

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

    if(forking == TRUE){
      net <- parallel::mclapply(net, function(x) add_setdiff(x), mc.cores = ncpus)
    } else {
      net <- lapply(net, function(x) add_setdiff(x))
    }

    vertices <- vertices[order(vertices)]

    N = max(unlist(lapply(net, function(x) network::network.size(x))))
  } else {
    vertices <- unique(unlist(lapply(net, function(x) network::get.vertex.attribute(x, 'vertex.names'))))
    vertices <- vertices[order(vertices)]
  }

  if(forking == TRUE){
    YT <- parallel::mclapply(net, function(x) network::as.matrix.network(x), mc.cores = ncpus)
  } else{
    YT <- lapply(net, function(x) network::as.matrix.network(x))
  }

  ### find each ego-network; use K steps out from each node
  # use Y+t(Y) to jump backwards across links at the second step out

  if(directed == FALSE){
    YT2 <- lapply(YT, function(x) x+t(x))
    recode <- function(x){
      x[x == 2] <- 1
      x <- x[order(colnames(x)),order(colnames(x))]
      return(x)
    }
    if(forking == TRUE){
      YT2 <- parallel::mclapply(YT2, function(x) recode(x), mc.cores = ncpus)
    } else {
      YT2 <- lapply(YT2, function(x) recode(x))
    }
  }

  if(directed == TRUE){
    YT2 <- YT
    recode <- function(x){
      x[x == 2] <- 1
      x <- x[order(colnames(x)),order(colnames(x))]
      return(x)
    }
    if(forking == TRUE){
      YT2 <- parallel::mclapply(YT2, function(x) recode(x), mc.cores = ncpus)
    } else{
      YT2 <- lapply(YT2, function(x) recode(x))
    }
  }

  neighborhood_extract <- function(x){
    net_x <- sna::gapply(x,c(1,2),1:N,"*",1,distance=core_size)
  }


  # reorder matrix by name, preserving order
  if(forking == TRUE){
    net2 <- parallel::mclapply(YT2, function(x) network::network(x, directed = directed), mc.cores = ncpus) #network object based upon the network matrix y which takes y and transforms it by	causing nodes to "jump backwards across links at the second step"
    xt<-parallel::mclapply(net2, function(x) neighborhood_extract(x), mc.cores = ncpus)
  } else {
    net2 <- lapply(YT2, function(x) network::network(x, directed = directed)) #network object based upon the network matrix y which takes y and transforms it by	causing nodes to "jump backwards across links at the second step"
    xt<-lapply(net2, function(x) neighborhood_extract(x))
  }

  #The following for loop says for every node. apply a function over the list of vectors. I think this is just a way to change the matrix into a more readable x,y form
  # reduce_adjacency function
  reduce_adjacency <- function(i){
    out <- as.matrix(Y[c(i,time_period[[i]]),c(i,time_period[[i]])])
    return(out)
  }

  # for each time period ti in time_steps
  for(ti in 1:time_steps){
    # extract the neighborhoods of each node for that time period
    time_period <- xt[[ti]]
    # create a new empty list x which will contain all the new ego-networks
    x <- list()
    # remove from TY2 the adjacency matrix associated with time period ti
    Y <- YT2[[ti]]
    # loop over each node from time_period's list of node-based neighborhoods and reduce the
    # broader adjacency matrix Y to the ego-network's adjacency matrix and save that in n
    if(forking == TRUE){
      x<-parallel::mclapply(seq_along(time_period), reduce_adjacency, mc.cores = ncpus)
      # make all the adjacency matrices into network objects
      x<-parallel::mclapply(x, function(x) network::as.network.matrix(x, directed=directed), mc.cores = ncpus)
    } else {
      x<-lapply(seq_along(time_period), reduce_adjacency)
      x<-lapply(x, function(x) network::as.network.matrix(x, directed=directed))
    }
    # overwrite ti's spot in the list of neighborhood indices with al the period's new ego-networks
    xt[[ti]] <- x
  }

  # get this in xt[[ego]][[time-observation]] format,
  rearrange_format <- function(i){
    time_list <- list()
    for(ti in 1:time_steps){
      nit <- xt[[ti]][[i]]
      time_list[[ti]] <- nit
    }
    return(time_list)
  }

  if(forking == TRUE){
    xt <- parallel::mclapply(as.list(1:N), rearrange_format, mc.cores =2)
  } else {
    xt <- lapply(as.list(1:N), rearrange_format)
  }

  # only look at egos bigger than MINSIZE.  Here I think the MINSIZE is however many connections they have
  keep_mat <- matrix(ncol = length(net), nrow = N)
  colnames(keep_mat) <- paste0("t", 1:time_steps)
  rownames(keep_mat) <- vertices[order(vertices)]

  # create function that evaluates the size of each ego-network and determine if they're large enough to be included
  for(i in 1:length(xt)){
    ego <- xt[[i]]
    keep <- lapply(ego, network::network.size)>=min_size
    keep_mat[i,] <- keep
  }

  # Unname the list
  red_net_list <- list()
  if(length(names(net)) > 0){
    net <- unname(net)
  }

  # Double check to make sure I'm removing those nodes whose ego-networks do not achieve the minimum size
  for(t in 1:time_steps){
    # extract time slice of the network
    time_slice <- net[[t]]
    # convert to adjacency
    red_net <- network::as.sociomatrix(time_slice)
    # order the adjacency using R's ordering
    red_net <- red_net[order(colnames(red_net)),order(colnames(red_net))]

    # If the reduced sociomatrix is null or has no rows or zero rows, initalize an empty matrix,
    # otherise reduce the adjacency by the indices from the keep matrix
    if(is.null(nrow(red_net[keep_mat[,t],keep_mat[,t]])) || nrow(red_net[keep_mat[,t],keep_mat[,t]]) == 0){
      red_net <- NA
    } else {
      red_net <- network::delete.vertices(x = time_slice, vid = which(get.vertex.attribute(time_slice, 'vertex.names') %in% get.vertex.attribute(time_slice, 'vertex.names')[keep_mat[,t]==FALSE]))
    }

    # If the reduced network is of size zero
    #if(network::network.size(red_net) != 0){
    #  for(att in network::list.vertex.attributes(time_slice)){
    #
    #   vals <- network::get.vertex.attribute(time_slice,att)
    #    names(vals) <- network::get.vertex.attribute(time_slice, 'vertex.names')

    #        to_match <- network::get.vertex.attribute(red_net, 'vertex.names')
    #       red_vals <- vals[to_match]

    #        network::set.vertex.attribute(red_net,att,red_vals)
    #      }
    #      if(edge_covariates == TRUE){
    #        for(att_e in setdiff(network::list.network.attributes(time_slice), c("bipartite", "directed", "hyper", "loops", "mnext", "multiple", "n"))){
    #          # setting edge attributes harder given how they're stored
    #          if(att_e != "na"){
    #            adj <- as.matrix(time_slice, matrix.type = "adjacency")
    #            adj <- adj[order(as.integer(colnames(adj))),order(as.integer(colnames(adj)))]
    #
    #            net_att <- network::get.network.attribute(time_slice, att_e)
    #
    #            colnames(net_att) <- get.vertex.attribute(orig_nets[[t]], 'vertex.names')
    #            rownames(net_att) <- get.vertex.attribute(orig_nets[[t]], 'vertex.names')

    #            adj_red <- as.matrix(red_net, matrix.type = "adjacency")
    #            vertices_red <- rownames(adj_red)
    #            net_red <- net_att[vertices_red, vertices_red]
    #
    #            network::set.network.attribute(red_net,att_e,net_red)
    #
    #         }
    #        }
    #    }
    #    }
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

  orig_xt <- xt

  # Now i'm going to populate each of the new ego-networks with attributes from the broader networks
  populate_attributes <- function(i){
    new_nets <- orig_xt[[i]]
    for(t in 1:length(new_nets)){
      if(length(new_nets) > 1){
        vertex_ids <- network::get.vertex.attribute(new_nets[[t]], 'vertex.names')
        time_indices <- which(keep_mat[i,] == TRUE)
        index <- time_indices[t]
        indices <- which(network::get.vertex.attribute(orig_nets[[index]], 'vertex.names') %in% vertex_ids)
        for(att in network::list.vertex.attributes(orig_nets[[index]])){
          network::set.vertex.attribute(new_nets[[t]], att, network::get.vertex.attribute(orig_nets[[index]], att)[indices])
        }
        if(edge_covariates == TRUE){
          for(att_e in setdiff(network::list.network.attributes(orig_nets[[index]]), c("bipartite", "directed", "hyper", "loops", "mnext", "multiple", "n"))){
            if(att_e != "na"){

              adj <- as.matrix(orig_nets[[index]], matrix.type = "adjacency")
              net_att <- network::get.network.attribute(orig_nets[[index]], att_e)
              colnames(net_att) <- colnames(adj)[1:nrow(net_att)]
              rownames(net_att) <- rownames(adj)[1:nrow(net_att)]

              adj_red <- as.matrix(new_nets[[t]], matrix.type = "adjacency")
              vertices_red <- stats::na.omit(rownames(adj_red))
              net_red <- net_att[vertices_red, vertices_red]

              network::set.network.attribute(new_nets[[t]],att_e,net_red)
            }
          }
        }
      }
      if(length(new_nets) == 1 & t == 1){
        vertex_ids <- network::get.vertex.attribute(new_nets[[1]], 'vertex.names')
        time_indices <- which(keep_mat[i,] == TRUE)
        index <- time_indices[1]
        indices <- which(network::get.vertex.attribute(orig_nets[[index]], 'vertex.names') %in% vertex_ids) #
        for(att in network::list.vertex.attributes(orig_nets[[index]])){
          if(att != "na"){
            network::set.vertex.attribute(new_nets[[1]], att, network::get.vertex.attribute(orig_nets[[index]], att)[indices])
          }
        }
        if(edge_covariates == TRUE){
          for(att_e in setdiff(network::list.network.attributes(time_slice), c("bipartite", "directed", "hyper", "loops", "mnext", "multiple", "n"))){
            if(att_e != "na"){

              adj <- as.matrix(orig_nets[[index]], matrix.type = "adjacency")
              net_att <- network::get.network.attribute(orig_nets[[index]], att_e)
              colnames(net_att) <- colnames(adj)[1:nrow(net_att)]
              rownames(net_att) <- rownames(adj)[1:nrow(net_att)]

              adj_red <- as.matrix(new_nets[[1]], matrix.type = "adjacency")
              vertices_red <- stats::na.omit(rownames(adj_red))
              net_red <- net_att[vertices_red, vertices_red]

              network::set.network.attribute(new_nets[[1]],att_e,net_red)
            }
          }
        }
      }
      if(length(new_nets) == 0){
        new_nets <- NA
      }
    }
    return(new_nets)
  }

  if(forking == TRUE){
    xt <- parallel::mclapply(seq_along(orig_xt), populate_attributes, mc.cores = ncpus)
  } else {
    xt <- lapply(seq_along(orig_xt), populate_attributes)
  }

  null_networks <- which(is.na(xt))
  dropped <- vertices[null_networks]
  if(length(dropped) > 0){
    cat("Longitudinally observed ego-networks removed because no time periods achieve the minimum size necessary to be included and to allow for an identifiable model.  If necessary, adjust min_size parameter.", "\n")
    print(paste(dropped, "removed from ego-TERGM analysis due to no identifiable networks"))
  }

  remaining_vertices <- vertices[-c(null_networks)]

  xt <- xt[!is.na(xt)]

  x <- xt
  rm(xt)
  cat("Data formatting complete.", "\n")

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

  # Variation on Leifeld's btergm function to overcome errors from separation.
  createBtergm_local <-function(coef, boot, R, nobs, time.steps, formula, formula2,
                                response, effects, weights, auto.adjust, offset, directed,
                                bipartite, nvertices, data){
    new("btergm", coef = coef, boot = boot, R = R, nobs = nobs,
        time.steps = time.steps, formula = formula, formula2 = formula2,
        response = response, effects = effects, weights = weights,
        auto.adjust = auto.adjust, offset = offset, directed = directed,
        bipartite = bipartite, nvertices = nvertices, data = data)
  }

  tergmprepare_local <- function (formula, offset = TRUE, blockdiag = FALSE, verbose = TRUE){
    l <- list()
    l$lhs.original <- deparse(formula[[2]])
    l$networks <- eval(parse(text = deparse(formula[[2]])), envir = environment(formula))
    if (class(l$networks) == "list" || class(l$networks) == "network.list") {
    }
    else {
      l$networks <- list(l$networks)
    }
    for (i in 1:length(l$networks)) {
      if (!class(l$networks[[i]]) %in% c("network", "matrix",
                                         "list")) {
        tryCatch({
          l$networks[[i]] <- as.matrix(l$networks[[i]])
        }, error = function(cond) {
          stop(paste("Object", i, "could not be converted to a matrix."))
        })
      }
    }
    l$num.vertices <- max(sapply(l$networks, function(x) network::get.network.attribute(network::network(x),
                                                                                        "n")))
    if (network::is.network(l$networks[[1]])) {
      l$directed <- network::is.directed(l$networks[[1]])
      l$bipartite <- network::is.bipartite(l$networks[[1]])
    }
    else {
      if (xergm.common::is.mat.directed(as.matrix(l$networks[[1]]))) {
        l$directed <- TRUE
      }
      else {
        l$directed <- FALSE
      }
      if (xergm.common::is.mat.onemode(as.matrix(l$networks[[1]]))) {
        l$bipartite <- FALSE
      }
      else {
        l$bipartite <- TRUE
      }
    }
    l$form <- update.formula(formula, networks[[i]] ~ .)
    l$time.steps <- length(l$networks)
    tilde <- deparse(l$form[[1]])
    lhs <- deparse(l$form[[2]])
    rhs <- paste(deparse(l$form[[3]]), collapse = "")
    rhs <- gsub("\\s+", " ", rhs)
    rhsterms <- strsplit(rhs, "\\s*\\+\\s*")[[1]]
    if (length(rhsterms) > 1) {
      for (i in length(rhsterms):2) {
        leftbracketmatches <- gregexpr("\\(", rhsterms[i])[[1]]
        leftbracketmatches <- leftbracketmatches[leftbracketmatches !=
                                                   -1]
        leftbracketmatches <- length(leftbracketmatches)
        rightbracketmatches <- gregexpr("\\)", rhsterms[i])[[1]]
        rightbracketmatches <- rightbracketmatches[rightbracketmatches !=
                                                     -1]
        rightbracketmatches <- length(rightbracketmatches)
        if (leftbracketmatches != rightbracketmatches) {
          rhsterms[i - 1] <- paste(rhsterms[i - 1], rhsterms[i],
                                   sep = " + ")
          rhsterms <- rhsterms[-i]
        }
      }
    }
    l$rhs.terms <- rhsterms
    rhs.operators <- rep("+", length(l$rhs.terms) - 1)
    covnames <- character()
    for (k in 1:length(l$rhs.terms)) {
      if (grepl("((edge)|(dyad))cov", l$rhs.terms[k])) {
        if (grepl(",\\s*?((attr)|\\\")", l$rhs.terms[k])) {
          s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)"
        }
        else {
          s <- "((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)"
        }
        x1 <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
        x2 <- sub(s, "\\5", l$rhs.terms[k], perl = TRUE)
        if (grepl("\\[.*\\]", x2)) {
          stop(paste0("Covariate names are not allowed to have indices: ",
                      x2, ". Please prepare a list object before estimation."))
        }
        if (grepl("^\"", x2))
          next
        x3 <- sub(s, "\\6", l$rhs.terms[k], perl = TRUE)
        x.current <- eval(parse(text = x2), envir = environment(formula))
        type <- class(x.current)
        l$covnames <- c(l$covnames, x2)
        l[[x2]] <- x.current
        if (grepl("\\[i\\]+$", x2)) {
          stop(paste0("Error in the following model term: ",
                      l$rhs.terms[k], ". The index 'i' is used internally by btergm. Please use a ",
                      "different index, for example 'j'."))
        }
        if (grepl("[^\\]]\\]$", x2)) {
          l$rhs.terms[k] <- paste0(x1, x2, x3)
          if (type %in% c("matrix", "network", "dgCMatrix",
                          "dgTMatrix", "dsCMatrix", "dsTMatrix", "dgeMatrix")) {
            x.current <- list(x.current)
            l[[x2]] <- x.current
          }
          if (length(x.current) != l$time.steps) {
            stop(paste(x2, "has", length(x.current), "elements, but there are",
                       l$time.steps, "networks to be modeled."))
          }
          if (blockdiag == TRUE) {
          }
          else {
            x2 <- paste0(x2, "[[i]]")
          }
        }
        else if (type %in% c("matrix", "network", "dgCMatrix",
                             "dgTMatrix", "dsCMatrix", "dsTMatrix", "dgeMatrix")) {
          if (!type %in% c("matrix", "network")) {
            x.current <- as.matrix(x.current)
          }
          l[[x2]] <- list()
          for (i in 1:l$time.steps) {
            l[[x2]][[i]] <- x.current
          }
          if (blockdiag == TRUE) {
          }
          else {
            x2 <- paste0(x2, "[[i]]")
          }
          l$rhs.terms[k] <- paste(x1, x2, x3, sep = "")
        }
        else if (type == "list" || type == "network.list") {
          if (length(x.current) != l$time.steps) {
            stop(paste(x2, "has", length(get(x2)), "elements, but there are",
                       l$time.steps, "networks to be modeled."))
          }
          if (blockdiag == TRUE) {
          }
          else {
            x2 <- paste0(x2, "[[i]]")
          }
          l$rhs.terms[k] <- paste0(x1, x2, x3)
        }
        else {
          tryCatch({
            l[[x2]] <- list(rep(as.matrix(x.current)),
                            l$time.steps)
          }, error = function(cond) {
            stop(paste0("Object '", x2, "' could not be converted to a matrix."))
          })
        }
      }
      else if (grepl("memory", l$rhs.terms[k])) {
        s <- "(?:memory\\((?:.*type\\s*=\\s*)?(?:\"|'))(\\w+)(?:(\"|').*\\))"
        if (grepl(s, l$rhs.terms[k]) == FALSE) {
          type <- "stability"
        }
        else {
          type <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
        }
        s <- "(?:memory\\(.*lag\\s*=\\s*)(\\d+)(?:.*\\))"
        if (grepl(s, l$rhs.terms[k]) == FALSE) {
          lag <- 1
        }
        else {
          lag <- as.integer(sub(s, "\\1", l$rhs.terms[k],
                                perl = TRUE))
        }
        if (lag > length(l$networks) - 1) {
          stop("The 'lag' argument in the 'memory' term is too large.")
        }
        mem <- l$networks[-(length(l$networks):(length(l$networks) -
                                                  lag + 1))]
        mem <- lapply(mem, as.matrix)
        memory <- list()
        for (i in 1:length(mem)) {
          if (type == "autoregression") {
            memory[[i]] <- mem[[i]]
          }
          else if (type == "stability") {
            mem[[i]][mem[[i]] == 0] <- -1
            memory[[i]] <- mem[[i]]
          }
          else if (type == "innovation") {
            memory[[i]] <- mem[[i]]
            memory[[i]][mem[[i]] == 0] <- 1
            memory[[i]][mem[[i]] == 1] <- 0
          }
          else if (type == "loss") {
            memory[[i]] <- mem[[i]]
            memory[[i]][mem[[i]] == 0] <- 0
            memory[[i]][mem[[i]] == 1] <- -1
          }
          else {
            stop("'type' argument in the 'memory' term not recognized.")
          }
        }
        rm(mem)
        l[["memory"]] <- memory
        if (blockdiag == TRUE) {
          l$rhs.terms[k] <- "edgecov(memory)"
        }
        else {
          l$rhs.terms[k] <- "edgecov(memory[[i]])"
        }
        l$covnames <- c(l$covnames, "memory")
      }
      else if (grepl("delrecip", l$rhs.terms[k])) {
        s <- "(?:delrecip\\((?:.*mutuality\\s*=\\s*)?)((TRUE)|(FALSE)|T|F)(?:.*\\))"
        if (grepl(s, l$rhs.terms[k]) == FALSE) {
          mutuality <- FALSE
        }
        else {
          mutuality <- as.logical(sub(s, "\\1", l$rhs.terms[k],
                                      perl = TRUE))
        }
        s <- "(?:delrecip\\(.*lag\\s*=\\s*)(\\d+)(?:.*\\))"
        if (grepl(s, l$rhs.terms[k]) == FALSE) {
          lag <- 1
        }
        else {
          lag <- as.integer(sub(s, "\\1", l$rhs.terms[k],
                                perl = TRUE))
        }
        if (lag > length(l$networks) - 1) {
          stop("The 'lag' argument in the 'delrecip' term is too large.")
        }
        dlr <- l$networks[-(length(l$networks):(length(l$networks) -
                                                  lag + 1))]
        dlr <- lapply(dlr, function(x) t(as.matrix(x)))
        delrecip <- list()
        for (i in 1:length(dlr)) {
          delrecip[[i]] <- dlr[[i]]
          if (mutuality == TRUE) {
            delrecip[[i]][dlr[[i]] == 0] <- -1
          }
        }
        rm(dlr)
        l[["delrecip"]] <- delrecip
        if (blockdiag == TRUE) {
          l$rhs.terms[k] <- "edgecov(delrecip)"
        }
        else {
          l$rhs.terms[k] <- "edgecov(delrecip[[i]])"
        }
        l$covnames <- c(l$covnames, "delrecip")
      }
      else if (grepl("timecov", l$rhs.terms[k])) {
        s <- "(?:timecov\\((?:.*x\\s*=\\s*)?)(\\w+)(?:.*\\))"
        if (sub(s, "\\1", l$rhs.terms[k], perl = TRUE) %in%
            c("minimum", "maximum", "transform", "min", "max",
              "trans")) {
          s <- "(?:timecov\\(?:.*x\\s*=\\s*)(\\w+)(?:.*\\))"
        }
        countprevtc <- 1
        if (k > 1) {
          for (i in (k - 1):1) {
            if (grepl("timecov", l$rhs.terms[i])) {
              countprevtc <- countprevtc + 1
            }
          }
        }
        if (countprevtc > 0) {
          countprevtc <- as.character(countprevtc)
        }
        else {
          countprevtc <- ""
        }
        if (grepl(s, l$rhs.terms[k]) == FALSE) {
          x <- NULL
          label <- paste0("timecov", countprevtc)
        }
        else {
          x <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
          label <- paste0("timecov", countprevtc, ".",
                          x)
        }
        s <- "(?:timecov\\(.*minimum\\s*=\\s*)(\\d+)(?:.*\\))"
        if (grepl(s, l$rhs.terms[k]) == FALSE) {
          minimum <- 1
        }
        else {
          minimum <- as.integer(sub(s, "\\1", l$rhs.terms[k],
                                    perl = TRUE))
        }
        s <- "(?:timecov\\(.*maximum\\s*=\\s*)(\\d+)(?:.*\\))"
        if (grepl(s, l$rhs.terms[k]) == FALSE) {
          maximum <- l$time.steps
        }
        else {
          maximum <- as.integer(sub(s, "\\1", l$rhs.terms[k],
                                    perl = TRUE))
        }
        s <- "(?:timecov\\(.*transform\\s*=\\s*)(.+?)(?:(?:,|\\)$)]*.*)"
        if (grepl(s, l$rhs.terms[k]) == FALSE) {
          transform <- function(t) t
        }
        else {
          transform <- eval(parse(text = sub(s, "\\1",
                                             l$rhs.terms[k], perl = TRUE)))
        }
        if (is.null(x)) {
          covariate <- l[["networks"]]
          onlytime <- TRUE
        }
        else {
          onlytime <- FALSE
          covariate <- get(x)
        }

        timecov <- function(covariate, minimum = 1, maximum = length(covariate),
                            transform = function(t) 1 + (0 * t) + (0 * t^2), onlytime = FALSE)
        {
          if (class(covariate) != "list") {
            stop("'covariate' must be a list of matrices or network objects.")
          }
          for (i in 1:length(covariate)) {
            if (is.network(covariate[[i]])) {
              covariate[[i]] <- as.matrix(covariate[[i]])
            }
            else if (!is.matrix(covariate[[i]])) {
              stop("'covariate' must be a list of matrices or network objects.")
            }
          }
          if (is.null(minimum) || is.null(maximum) || !is.numeric(minimum) ||
              !is.numeric(maximum) || length(minimum) > 1 || length(maximum) >
              1) {
            stop("'minimum' and 'maximum' must be single numeric values.")
          }
          if (is.null(transform)) {
            transform <- function(t) 1 + (0 * t) + (0 * t^2)
          }
          else if (!is.function(transform)) {
            stop("'transform' must be a function.")
          }
          l <- 1:length(covariate)
          values <- transform(l)
          if (is.null(values) || any(is.null(values)) || any(!is.finite(values)) ||
              any(is.na(values)) || any(!is.numeric(values))) {
            stop("The 'transform' function produces non-numeric values.")
          }
          values <- values * (l >= minimum) * (l <= maximum)
          timecov <- list()
          for (i in 1:length(l)) {
            if (onlytime == FALSE) {
              timecov[[i]] <- covariate[[i]] * matrix(values[i],
                                                      nrow = nrow(covariate[[i]]), ncol = ncol(covariate[[i]]))
            }
            else {
              timecov[[i]] <- matrix(values[i], nrow = nrow(covariate[[i]]),
                                     ncol = ncol(covariate[[i]]))
            }
          }
          for (i in 1:length(timecov)) {
            rownames(timecov[[i]]) <- rownames(covariate[[i]])
            colnames(timecov[[i]]) <- colnames(covariate[[i]])
          }
          return(timecov)
        }
        tc <- timecov(covariate = covariate, minimum = minimum,
                      maximum = maximum, transform = transform, onlytime = onlytime)
        l[[label]] <- tc
        labelsuffix <- sub(s, "\\1", l$rhs.terms[k], perl = TRUE)
        labelsuffix <- if (blockdiag == TRUE) {
          l$rhs.terms[k] <- paste0("edgecov(", label, ")")
        }
        else {
          l$rhs.terms[k] <- paste0("edgecov(", label, "[[i]])")
        }
        l$covnames <- c(l$covnames, label)
        if (verbose == TRUE) {
          timecovreporting <- matrix(sapply(tc, function(x) mean(x[1,
                                                                   2])), nrow = 1)
          colnames(timecovreporting) <- paste0("t=", 1:length(l$networks))
          rownames(timecovreporting) <- ""
          message("Mean transformed timecov values:")
          print(timecovreporting)
        }
      }
    }
    l$covnames <- c("networks", l$covnames)
    lengths <- sapply(l$covnames, function(cn) length(l[[cn]]))
    mn <- max(lengths)
    if (length(table(lengths)) > 1) {
      mn <- min(lengths)
      l$time.steps <- mn
      for (i in 1:length(l$covnames)) {
        cn <- l$covnames[[i]]
        ll <- l[[cn]]
        difference <- length(ll) - mn
        if (difference > 0) {
          l[[cn]] <- ll[(difference + 1):length(ll)]
        }
      }
    }
    t.end <- max(lengths)
    t.start <- t.end - mn + 1
    if (verbose == TRUE) {
      if (length(l$covnames) > 1) {
        dimensions <- lapply(lapply(l$covnames, function(x) l[[x]]),
                             function(y) sapply(y, function(z) dim(as.matrix(z))))
        rownames(dimensions[[1]]) <- paste(l$lhs.original,
                                           c("(row)", "(col)"))
        for (i in 2:length(dimensions)) {
          rownames(dimensions[[i]]) <- c(paste(l$covnames[i],
                                               "(row)"), paste(l$covnames[i], "(col)"))
        }
        dimensions <- do.call(rbind, dimensions)
        colnames(dimensions) <- paste0("t=", t.start:t.end)
        message("\nInitial dimensions of the network and covariates:")
        print(dimensions)
      }
      else {
        message("\nNo covariates provided.")
      }
    }
    l$auto.adjust <- FALSE
    if (length(l$covnames) > 1) {
      nr <- lapply(lapply(l$covnames, function(x) l[[x]]),
                   function(y) sapply(y, function(z) nrow(as.matrix(z))))
      nr <- do.call(rbind, nr)
      nc <- lapply(lapply(l$covnames, function(x) l[[x]]),
                   function(y) sapply(y, function(z) ncol(as.matrix(z))))
      nc <- do.call(rbind, nc)
      for (i in 1:ncol(nr)) {
        if (length(unique(nr[, i])) > 1) {
          l$auto.adjust <- TRUE
        }
      }
      for (i in 1:ncol(nc)) {
        if (length(unique(nc[, i])) > 1) {
          l$auto.adjust <- TRUE
        }
      }
      if (verbose == TRUE && l$auto.adjust == TRUE) {
        message(paste("\nDimensions differ across networks within time steps."))
      }
      if (l$auto.adjust == TRUE) {
        for (i in 1:length(l$covnames)) {
          for (t in 1:l$time.steps) {
            if (is.null(rownames(as.matrix(l[[l$covnames[i]]][[t]]))) ||
                is.null(colnames(as.matrix(l[[l$covnames[i]]][[t]])))) {
              stop(paste0("The dimensions of the covariates differ, but ",
                          "covariate '", l$covnames[i], " does not have node labels at t = ",
                          t, ". Automatic adjustment of dimensions is therefore not ",
                          "possible."))
            }
          }
        }
      }
      if (l$auto.adjust == FALSE) {
        for (t in 1:l$time.steps) {
          rlabels.i <- list()
          clabels.i <- list()
          for (i in 1:length(l$covnames)) {
            rlabels.i[[i]] <- rownames(as.matrix(l[[l$covnames[i]]][[t]]))
            clabels.i[[i]] <- colnames(as.matrix(l[[l$covnames[i]]][[t]]))
          }
          rlabels.i <- do.call(rbind, rlabels.i)
          clabels.i <- do.call(rbind, clabels.i)
          flag <- FALSE
          if (!is.null(rlabels.i)) {
            for (j in 1:ncol(rlabels.i)) {
              if (length(unique(rlabels.i[, j])) > 1) {
                l$auto.adjust <- TRUE
                flag <- TRUE
                break
              }
            }
          }
          if (!is.null(clabels.i)) {
            for (j in 1:ncol(clabels.i)) {
              if (length(unique(clabels.i[, j])) > 1) {
                l$auto.adjust <- TRUE
                flag <- TRUE
                break
              }
            }
          }
        }
        if (verbose == TRUE && flag == TRUE) {
          message(paste("\nSame dimensions but different labels across",
                        "networks within time steps."))
        }
      }
    }
    if (verbose == TRUE && l$auto.adjust == TRUE) {
      message("Trying to auto-adjust the dimensions of the networks. ",
              "If this fails, provide conformable matrices or network objects.")
    }
    else if (verbose == TRUE) {
      message("\nAll networks are conformable.")
    }
    structzero.df <- data.frame(label = character(), time = integer(),
                                object = character(), where = character())
    if (length(l$covnames) > 0 && l$auto.adjust == TRUE) {
      for (i in 1:l$time.steps) {
        for (j in 1:length(l$covnames)) {
          for (k in 1:length(l$covnames)) {
            if (j != k) {
              nw.j <- l[[l$covnames[j]]][[i]]
              rn.j <- rownames(as.matrix(nw.j))
              cn.j <- colnames(as.matrix(nw.j))
              nr.j <- nrow(as.matrix(nw.j))
              nc.j <- ncol(as.matrix(nw.j))
              nw.k <- l[[l$covnames[k]]][[i]]
              rn.k <- rownames(as.matrix(nw.k))
              cn.k <- colnames(as.matrix(nw.k))
              nr.k <- nrow(as.matrix(nw.k))
              nc.k <- ncol(as.matrix(nw.k))
              if (is.null(rn.j) || is.null(cn.j)) {
                stop(paste0("Missing row or column labels in object '",
                            l$covnames[j], "'. Provide row and column ",
                            "labels for all networks and covariates."))
              }
              else if (is.null(rn.k) || is.null(cn.k)) {
                stop(paste0("Missing row or column labels in object '",
                            l$covnames[k], "'. Provide row and column ",
                            "labels for all networks and covariates."))
              }
              else {
                if (is.null(rn.j) && !is.null(rn.k) &&
                    nr.j == nr.k) {
                  if (class(nw.j) == "network") {
                    network::set.vertex.attribute(nw.j,
                                                  "vertex.names", rn.k)
                  }
                  else {
                    rownames(nw.j) <- rn.k
                  }
                }
                else if (is.null(rn.k) && !is.null(rn.j) &&
                         nr.j == nr.k) {
                  if (class(nw.k) == "network") {
                    network::set.vertex.attribute(nw.k,
                                                  "vertex.names", rn.j)
                  }
                  else {
                    rownames(nw.k) <- rn.j
                  }
                }
                else if ((is.null(rn.k) || is.null(rn.j)) &&
                         nr.j != nr.k) {
                  stop(paste0("Object '", l$covnames[j],
                              "' is incompatible with object '",
                              l$covnames[k], "' at t = ", i, "."))
                }
                nw.j.labels <- adjust(nw.j, nw.k, remove = FALSE,
                                      value = 1, returnlabels = TRUE)
                nw.j <- adjust(nw.j, nw.k, remove = FALSE,
                               value = 1)
                l[[l$covnames[j]]][[i]] <- nw.j
                ro <- nw.j.labels$added.row
                co <- nw.j.labels$added.col
                if (length(ro) > 0) {
                  ro <- data.frame(label = ro, time = rep(i,
                                                          length(ro)), object = rep(l$covnames[j],
                                                                                    length(ro)), where = rep("row", length(ro)))
                  structzero.df <- rbind(structzero.df,
                                         ro)
                }
                if (length(co) > 0) {
                  co <- data.frame(label = co, time = rep(i,
                                                          length(co)), object = rep(l$covnames[j],
                                                                                    length(co)), where = rep("col", length(co)))
                  structzero.df <- rbind(structzero.df,
                                         co)
                }
                nw.k.labels <- adjust(nw.k, nw.j, remove = FALSE,
                                      value = 1, returnlabels = TRUE)
                nw.k <- adjust(nw.k, nw.j, remove = FALSE,
                               value = 1)
                l[[l$covnames[k]]][[i]] <- nw.k
                ro <- nw.k.labels$added.row
                co <- nw.k.labels$added.col
                if (length(ro) > 0) {
                  ro <- data.frame(label = ro, time = rep(i,
                                                          length(ro)), object = rep(l$covnames[j],
                                                                                    length(ro)), where = rep("row", length(ro)))
                  structzero.df <- rbind(structzero.df,
                                         ro)
                }
                if (length(co) > 0) {
                  co <- data.frame(label = co, time = rep(i,
                                                          length(co)), object = rep(l$covnames[j],
                                                                                    length(co)), where = rep("col", length(co)))
                  structzero.df <- rbind(structzero.df,
                                         co)
                }
              }
            }
          }
        }
      }
    }
    nr.net <- sapply(l$networks, function(x) nrow(as.matrix(x)))
    for (i in 1:length(l$covnames)) {
      nr <- sapply(l[[l$covnames[i]]], function(x) {
        nrow(as.matrix(x))
      })
      for (j in 1:l$time.steps) {
        if (nr[j] != nr.net[j]) {
          stop(paste0("Covariate object '", l$covnames[i],
                      "' does not have the same number of rows as the dependent ",
                      "network at time step ", j, "."))
        }
      }
    }
    nc.net <- sapply(l$networks, function(x) ncol(as.matrix(x)))
    for (i in 1:length(l$covnames)) {
      nc <- sapply(l[[l$covnames[i]]], function(x) {
        ncol(as.matrix(x))
      })
      for (j in 1:l$time.steps) {
        if (nc[j] != nc.net[j]) {
          stop(paste0("Covariate object '", l$covnames[i],
                      "' does not have the same number of columns as the dependent ",
                      "network at time step ", j, "."))
        }
      }
    }
    if (verbose == TRUE) {
      if (l$auto.adjust == TRUE) {
        sz.row <- unique(structzero.df[structzero.df$where ==
                                         "row", -3])
        szrownum <- numeric(length(l$networks))
        for (i in 1:length(l$networks)) {
          szrownum[i] <- nrow(sz.row[sz.row$time == i,
                                     ])
        }
        sz.col <- unique(structzero.df[structzero.df$where ==
                                         "col", -3])
        szcolnum <- numeric(length(l$networks))
        for (i in 1:length(l$networks)) {
          szcolnum[i] <- nrow(sz.col[sz.col$time == i,
                                     ])
        }
        totrow <- sapply(l$networks, function(x) nrow(as.matrix(x)))
        totcol <- sapply(l$networks, function(x) ncol(as.matrix(x)))
        if (offset == TRUE) {
          dimensions <- rbind(totrow, totcol, szrownum,
                              szcolnum, totrow - szrownum, totcol - szcolnum)
          rownames(dimensions) <- c("total number of rows",
                                    "total number of columns", "row-wise structural zeros",
                                    "column-wise structural zeros", "remaining rows",
                                    "remaining columns")
        }
        else {
          dimensions <- rbind(szrownum, szcolnum, totrow -
                                szrownum, totcol - szcolnum)
          rownames(dimensions) <- c("maximum deleted nodes (row)",
                                    "maximum deleted nodes (col)", "remaining rows",
                                    "remaining columns")
        }
        colnames(dimensions) <- paste0("t=", t.start:t.end)
        if (nrow(structzero.df) > 0) {
          if (offset == TRUE) {
            message("\nNodes affected completely by structural zeros:")
          }
          else {
            message("\nAbsent nodes:")
          }
          szcopy <- structzero.df
          szcopy$time <- szcopy$time - 1 + t.start
          print(unique(szcopy))
        }
        else {
          message("\nAll nodes are retained.")
        }
        message("\nNumber of nodes per time step after adjustment:")
        print(dimensions)
      }
    }
    l$nvertices <- sapply(l$networks, function(x) c(nrow(as.matrix(x)),
                                                    ncol(as.matrix(x))))
    rownames(l$nvertices) <- c("row", "col")
    colnames(l$nvertices) <- paste0("t=", t.start:t.end)
    l$offsmat <- list()
    for (i in 1:l$time.steps) {
      mat <- matrix(0, nrow = nrow(as.matrix(l$networks[[i]])),
                    ncol = ncol(as.matrix(l$networks[[i]])))
      rownames(mat) <- rownames(as.matrix(l$networks[[i]]))
      colnames(mat) <- colnames(as.matrix(l$networks[[i]]))
      l$offsmat[[i]] <- mat
    }
    if (nrow(structzero.df) > 0) {
      for (i in 1:nrow(structzero.df)) {
        if (structzero.df$where[i] == "row") {
          index <- which(rownames(l$offsmat[[structzero.df$time[i]]]) ==
                           structzero.df$label[i])
          l$offsmat[[structzero.df$time[i]]][index, ] <- 1
        }
        else {
          index <- which(colnames(l$offsmat[[structzero.df$time[i]]]) ==
                           structzero.df$label[i])
          l$offsmat[[structzero.df$time[i]]][, index] <- 1
        }
      }
    }
    if (offset == TRUE) {
      l$rhs.terms[length(l$rhs.terms) + 1] <- "offset(edgecov(offsmat[[i]]))"
      rhs.operators[length(rhs.operators) + 1] <- "+"
    }
    else {
      if (l$auto.adjust == TRUE) {
        l$offsmat <- suppressMessages(handleMissings(l$offsmat,
                                                     na = 1, method = "remove"))
        for (j in 1:length(l$covnames)) {
          l[[l$covnames[j]]] <- adjust(l[[l$covnames[j]]],
                                       l$offsmat)
        }
      }
    }
    if (verbose == TRUE && length(l$covnames) > 1) {
      dimensions <- lapply(lapply(l$covnames, function(x) l[[x]]),
                           function(y) sapply(y, function(z) dim(as.matrix(z))))
      rownames(dimensions[[1]]) <- paste(l$lhs.original, c("(row)",
                                                           "(col)"))
      for (i in 2:length(dimensions)) {
        rownames(dimensions[[i]]) <- c(paste(l$covnames[i],
                                             "(row)"), paste(l$covnames[i], "(col)"))
      }
      dimensions <- do.call(rbind, dimensions)
      colnames(dimensions) <- paste0("t=", t.start:t.end)
      message("\nDimensions of the network and covariates after adjustment:")
      print(dimensions)
    }
    rhs <- l$rhs.terms[1]
    if (length(rhs.operators) > 0) {
      for (i in 1:length(rhs.operators)) {
        rhs <- paste(rhs, rhs.operators[i], l$rhs.terms[i +
                                                          1])
      }
    }
    f <- paste(lhs, tilde, rhs)
    l$form <- as.formula(f, env = environment())
    if (blockdiag == TRUE) {
      if (l$bipartite == TRUE) {
        stop(paste("MCMC estimation is currently only supported for one-mode",
                   "networks. Use the btergm function instead."))
      }
      l$form <- update.formula(l$form, networks ~ .)
      l$form <- paste(deparse(l$form), collapse = "")
      l$form <- paste(l$form, "+ offset(edgecov(offsmat))")
      l$form <- as.formula(l$form, env = environment())
      if (length(l$covnames) > 1) {
        for (j in 2:length(l$covnames)) {
          l[[l$covnames[j]]] <- as.matrix(Matrix::bdiag(lapply(l[[l$covnames[j]]],
                                                               as.matrix)))
        }
      }
      l$offsmat <- as.matrix(Matrix::bdiag(l$offsmat))
      bdoffset <- lapply(l$networks, as.matrix)
      for (i in 1:length(bdoffset)) {
        bdoffset[[i]][, ] <- 1
      }
      bdoffset <- as.matrix((Matrix::bdiag(bdoffset) - 1) *
                              -1)
      l$offsmat <- l$offsmat + bdoffset
      rm(bdoffset)
      l$offsmat[l$offsmat > 0] <- 1
      if (class(l$networks[[1]]) == "network") {
        attrnames <- network::list.vertex.attributes(l$networks[[1]])
        attributes <- list()
        for (i in 1:length(l$networks)) {
          attrib <- list()
          for (j in 1:length(attrnames)) {
            attrib[[j]] <- network::get.vertex.attribute(l$networks[[i]],
                                                         attrnames[j])
          }
          attributes[[i]] <- attrib
          l$networks[[i]] <- as.matrix(l$networks[[i]])
        }
        l$networks <- network::network(as.matrix(Matrix::bdiag(l$networks)),
                                       directed = l$directed, bipartite = l$bipartite)
        for (i in 1:length(attrnames)) {
          attrib <- unlist(lapply(attributes, function(x) x[[i]]))
          network::set.vertex.attribute(l$networks, attrnames[i],
                                        attrib)
        }
      }
      else {
        l$networks <- network::network(as.matrix(Matrix::bdiag(l$networks)),
                                       directed = l$directed, bipartite = l$bipartite)
      }
      if (verbose == TRUE) {
        cat("\n")
      }
    }
    form3 <- paste(deparse(l$form[[3]]), collapse = "")
    form3 <- gsub("\\s+", " ", form3)
    l$form <- paste(deparse(l$form[[2]]), deparse(l$form[[1]]),
                    form3)
    return(l)
  }

  btergm_local <- function(formula, R = 500, offset = FALSE, returndata = FALSE,
                           parallel = c("no", "multicore", "snow"), ncpus = 1, cl = NULL,
                           verbose = TRUE, ...){
    l <- tergmprepare_local(formula = formula, offset = offset, verbose = verbose)
    for (i in 1:length(l$covnames)) {
      assign(l$covnames[i], l[[l$covnames[i]]])
    }
    assign("offsmat", l$offsmat)
    form <- as.formula(l$form)
    if (l$time.steps == 1) {
      warning(paste("The confidence intervals and standard errors are",
                    "meaningless because only one time step was provided."))
    }
    if (verbose == TRUE && returndata == FALSE) {
      if (parallel[1] == "no") {
        parallel.msg <- "on a single computing core"
      }
      else if (parallel[1] == "multicore") {
        parallel.msg <- paste("using multicore forking on",
                              ncpus, "cores")
      }
      else if (parallel[1] == "snow") {
        parallel.msg <- paste("using parallel processing on",
                              ncpus, "cores")
      }
      if (offset == TRUE) {
        offset.msg <- "with offset matrix and "
      }
      else {
        offset.msg <- "with "
      }
      message("\nStarting pseudolikelihood estimation ", offset.msg,
              R, " bootstrapping replications ", parallel.msg,
              "...")
    }
    else if(verbose == TRUE && returndata == TRUE){
      message("\nReturning data frame with change statistics.")
    }
    if(offset == TRUE){
      Y <- NULL
      X <- NULL
      W <- NULL
      O <- NULL
      for (i in 1:length(l$networks)) {
        nw <- ergm::ergm.getnetwork(form)
        model <- ergm::ergm.getmodel(form, nw, initialfit = TRUE)
        Clist <- ergm::ergm.Cprepare(nw, model)
        Clist.miss <- ergm::ergm.design(nw, model, verbose = FALSE)
        pl <- ergm::ergm.pl(Clist, Clist.miss, model, theta.offset = c(rep(FALSE,
                                                                           length(l$rhs.terms) - 1), TRUE), verbose = FALSE,
                            control = ergm::control.ergm(init = c(rep(NA,
                                                                      length(l$rhs.terms) - 1), 1)))
        Y <- c(Y, pl$zy[pl$foffset == 0])
        X <- rbind(X, cbind(data.frame(pl$xmat[pl$foffset ==
                                                 0, ], check.names = FALSE), i))
        W <- c(W, pl$wend[pl$foffset == 0])
        O <- c(O, pl$foffset[pl$foffset == 0])
      }
      term.names <- colnames(X)[-length(colnames(X))]
      term.names <- c(term.names, "time")
      colnames(X) <- term.names
    } else {
      Y <- NULL
      X <- NULL
      W <- NULL
      O <- NULL
      for (i in 1:length(l$networks)) {
        mpli <- ergm::ergmMPLE(form)
        Y <- c(Y, mpli$response)
        if (i > 1 && ncol(X) != ncol(mpli$predictor) + 1) {
          cn.x <- colnames(X)[-ncol(X)]
          cn.i <- colnames(mpli$predictor)
          names.x <- cn.x[which(!cn.x %in% cn.i)]
          names.i <- cn.i[which(!cn.i %in% cn.x)]
          if (length(names.x) > 0) {
            for (nm in 1:length(names.x)) {
              mpli$predictor <- cbind(mpli$predictor, rep(0,
                                                          nrow(mpli$predictor)))
              colnames(mpli$predictor)[ncol(mpli$predictor)] <- names.x[nm]
            }
          }
          if (length(names.i) > 0) {
            for (nm in 1:length(names.i)) {
              X <- cbind(X[, 1:(ncol(X) - 1)], rep(0, nrow(X)),
                         X[, ncol(X)])
              colnames(X)[ncol(X) - 1] <- names.i[nm]
            }
          }
        }
        X <- rbind(X, cbind(mpli$predictor, i))
        W <- c(W, mpli$weights)
      }
      term.names <- colnames(X)[-length(colnames(X))]
      term.names <- c(term.names, "time")
      X <- data.frame(X)
      colnames(X) <- term.names
    }
    unique.time.steps <- unique(X$time)
    x <- X[, -ncol(X)]
    x <- as.data.frame(x)
    if (returndata == TRUE) {
      return(cbind(Y, x))
    }
    xsparse <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
    if (ncol(xsparse) == 1) {
      stop("At least two model terms must be provided to estimate a TERGM.")
    }
    est <- speedglm::speedglm.wfit(y = Y, X = xsparse, weights = W, offset = O,
                                   family = binomial(link = logit), sparse = TRUE)
    startval <- coef(est)
    nobs <- est$n
    estimate <- function(unique.time.steps, bsi, Yi = Y, xsparsei = xsparse,
                         Wi = W, Oi = O, timei = X$time, startvali = startval) {
      indic <- unlist(lapply(bsi, function(x) which(timei ==
                                                      x)))
      tryCatch(expr = {
        return(coef(speedglm::speedglm.wfit(y = Yi[indic], X = xsparsei[indic,], weights = Wi[indic], offset = Oi[indic], family = binomial(link = logit),
                                            sparse = TRUE, start = startvali)))
      }, error = function(e) {
        return(coef(stats::glm.fit(y = Yi[indic], x = as.matrix(x)[indic, ], weights = Wi[indic], offset = Oi[indic], family = binomial(link = logit))))
        #  }, warning = function(w) {
        #   warning(w)
        #  }, finally = {
        #   coef(glm.fit(y = Yi[indic], x = as.matrix(x)[indic, ], weights = Wi[indic], offset = Oi[indic], family = binomial(link = logit)))
      })
    }
    coefs <- boot::boot(unique.time.steps, estimate, R = R, Yi = Y,
                        xsparsei = xsparse, Wi = W, Oi = O, timei = X$time, startvali = startval,
                        parallel = parallel, ncpus = ncpus, cl = cl)
    rm(X)
    if (ncol(coefs$t) == 1 && length(term.names) > 1 && coefs$t[1, 1] == "glm.fit: algorithm did not converge") {
      stop(paste("Algorithm did not converge. There might be a collinearity ",
                 "between predictors and/or dependent networks at one or more time",
                 "steps."))
    }
    colnames(coefs$t) <- term.names[1:(length(term.names) - 1)]
    names(startval) <- colnames(coefs$t)
    data <- list()
    for (i in 1:length(l$covnames)) {
      data[[l$covnames[i]]] <- l[[l$covnames[i]]]
    }
    data$offsmat <- l$offsmat
    btergm.object <- createBtergm_local(startval, coefs, R, nobs, l$time.steps,
                                        formula, l$form, Y, x, W, l$auto.adjust, offset, l$directed,
                                        l$bipartite, nvertices = l$nvertices, data)
    if (verbose == TRUE) {
      message("Done.")
    }
    return(btergm.object)
  }



  ########################################################################
  ### Initialization functions
  ########################################################################
  #Specify function in terms of ego.terms and G
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
        fit = btergm_local(formula = form, R=R, verbose = FALSE)@coef
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
  comparison_roles_out <- data.frame(Id = remaining_vertices,
                          Role = z)

  # comparison

  if(length(ego_tergm_fit$role_assignments$Role) == length(comparison_roles_out$Role)){
    cross_tab <- table(ego_tergm_fit$role_assignments$Role, comparison_roles_out$Role)
    rownames(cross_tab) <- paste0("Original Cluster: ", 1:nrow(cross_tab))
    colnames(cross_tab) <- paste0("Subset Cluster: ", 1:nrow(cross_tab))

  } else {
    cat("Note: Comparison sets incompatible.  There are nodes in the original model output that are not in the comparison model output,
        or there are nodes in the comparison model output that are not in the original model output.
        Cross-tabulation table unavailable, relying instead on comparison of proportions.
        Note that due to label switching, here groups are sorted according to value and compared according to similarity in proportions.", "\n")
    cross_tab <- rbind(unname(sort(prop.table(table(ego_tergm_fit$role_assignments$Role)))),
                       unname(sort(prop.table(table(comparison_roles_out$Role)))))

    rownames(cross_tab) <- c("Original Assignment Proportions", "Subset Assignment Proportions")

    colnames(cross_tab) <- paste0("Cluster ", 1:ncol(cross_tab))

  }

  cat("Done.", "\n")
  cat("Completed Time:", format(Sys.time(), "%a %b %d %X %Y"), "\n")
  return(list(model.fit = "egoTERGM Comparison",
              splitting_probability = splitting_probability,
              networks_sampled = which(sample_indices == 1),
              comparison_lambda = lambda,
              comparison_group.theta = group.theta,
              comparison_EE.BIC = EE.BIC,
              comparison_TS.BIC = TS.BIC,
              comparison_role_assignments = comparison_roles_out,
              comparison_ego_nets = x,
              comparison_table = cross_tab))

}
