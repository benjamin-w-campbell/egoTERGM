#' Estimation of ego-Temporal Exponential Random Graph Model (ego-TERGM) using Expectation Maximization (EM).
#'
#' This function estimates an ego-TERGM on a longitudinally observed network.  Currently the function does not support comparisons of whole networks.
#' @param net The longitudinally observed network that an ego-TERGM will be fit on.  Must be presented as a list of networks.  Any vertex attributes should be attached to networks.  Currently the function does not support comparisons of whole networks.
#' @param core_size The order of alters to include. The default value of one implies only looking at an ego's alters and the connections among them.
#' @param min_size  The minimum number of nodes an ego-network must achieve to be included.  Defaults to five.
#' @param roles The number of roles that should be fit.  Defaults to 3.
#' @param form The formula comprised of ERGM or TERGM terms used to distinguish between clusters assignments.  Specified as a vector of comma separated terms. No default.
#' @param add_drop Do nodes drop out of the network or enter it? If so, specify as the default TRUE.
#' @param directed Should the longitudinal network be treated as directed? If so, specify as the default TRUE.
#' @param edge_covariates Are edge covariates included in the form term? IF so, specify as TRUE.  No default.  These should be stored as network attributes.
#' @param seed The seed set to replicate analysis for pseudorandom number generator.
#' @param R The number of bootstrap replications that should be used for the estimation of a bootstrapped MPLE estimated TERGM for model initialization.  Defaults to 10.
#' @param parallel How the BTERGM should be computed, either with parallel processing ("multicore", "snow") or no parallel processing ("no").  Defaults to no.
#' @param ncpus The number of CPUs that should should be used for estimation, defaults to 1.
#' @param steps The number of default EM steps that should be taken, defaults to 50.
#' @param tol The difference in parameter estimates between EM iterations to determine if the algorithm has converged.  Defaults to 1e-6.
#' @return A list of model results, including lambda (the probability of assignments), group.theta (the roles by terms cluster centroids),
#'         EE.BIC (the Salter-Townshend and Murphy BIC cross-sectional BIC), TS.BIC (the Campbell BIC penalizing for time-steps),
#'        role_assignments (a data frame of the most likely assignments), reduced_networks (A list of the networks with excluded egos),
#'        ego_nets (a list of ego-networks), and ego_nets_used (N x T matrix of logicals here TRUE refers to ego-networks kept).
#' @keywords ego-TERGM
#' @references{
#'  Campbell, Benjamin W. (2018):
#'  Inferring Latent Roles in Longitudinal Networks.
#'  Forthcoming in \emph{Political Analysis}.
#'
#'  Leifeld, Philip, Skyler J. Cranmer and Bruce A. Desmarais (2017):
#'  Temporal Exponential Random Graph Models with btergm: Estimation and Bootstrap Confidence Intervals.
#'   \emph{Journal of Statistical Software} 83(6): 1-36. \url{http://dx.doi.org/10.18637/jss.v083.i06}
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
  cat("Start Time:", format(Sys.time(), "%a %b %d %X %Y"), "\n")

  cat("Data formatting started.", "\n")
  # tested prior to 9/3/17
  orig_nets <- net
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

    vertices <- vertices[order(as.integer(vertices))]


    N = max(unlist(lapply(net, function(x) network::network.size(x))))
  } else {
    vertices <- unique(unlist(lapply(net, function(x) network::get.vertex.attribute(x, 'vertex.names'))))
    vertices <- vertices[order(as.integer(vertices))]
  }


  YT <- lapply(net, function(x) network::as.matrix.network(x))

  ### find each ego-network; use K steps out from each node
  # use Y+t(Y) to jump backwards across links at the second step out

  if(directed == FALSE){
    YT2 <- lapply(YT, function(x) x+t(x))
    recode <- function(x){
      x[x == 2] <- 1
      x <- x[order(as.integer(colnames(x))),order(as.integer(colnames(x)))]
      return(x)
    }
    YT2 <- lapply(YT2, function(x) recode(x))

  }

  if(directed == TRUE){
    YT2 <- YT
    recode <- function(x){
      x[x == 2] <- 1
      x <- x[order(as.integer(colnames(x))),order(as.integer(colnames(x)))]
      return(x)
    }
    YT2 <- lapply(YT2, function(x) recode(x))

  }

  # reorder matrix by name, preserving order

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


  if (min_size<=1){
   stop("Minimum size must be greater than 1.  Please adjust the min_size argument to a value greater than 1.")
  }

  # only look at egos bigger than MINSIZE.  Here I think the MINSIZE is however many connections they have
  keep_mat <- matrix(ncol = length(net), nrow = N)
  colnames(keep_mat) <- paste0("t", 1:time_steps)
  rownames(keep_mat) <- vertices[order(as.integer(vertices))]
  for(i in 1:length(xt)){
    ego <- xt[[i]]
    keep <- lapply(ego, network::network.size)>=min_size
    keep_mat[i,] <- keep
  }


  ### good

  # all_net<-net

  # Unname the list
  red_net_list <- list()
  if(length(names(net)) > 0){
    net <- unname(net)
  }

  # Double check to make sure I'm removing those nodes whose ego-networks do not achieve the minimum size
  for(t in 1:length(net)){
    time_slice <- net[[t]]


    red_net <- network::as.sociomatrix(time_slice)
    red_net <- red_net[order(as.integer(colnames(red_net))),order(as.integer(colnames(red_net)))]

    if(is.null(nrow(red_net[keep_mat[,t],keep_mat[,t]])) || nrow(red_net[keep_mat[,t],keep_mat[,t]]) == 0){
      red_net <- network::network.initialize(n=0, directed=directed)
    } else {
      red_net <- network::as.network(red_net[keep_mat[,t],keep_mat[,t]], directed=directed)

    }

    if(network.size(red_net) == 0){
      red_net <- NA
    } else {
      for(att in network::list.vertex.attributes(time_slice)){

        vals <- network::get.vertex.attribute(time_slice,att)
        names(vals) <- network::get.vertex.attribute(time_slice, 'vertex.names')

        to_match <- network::get.vertex.attribute(red_net, 'vertex.names')
        red_vals <- vals[to_match]

        network::set.vertex.attribute(red_net,att,red_vals)
      }
      if(edge_covariates == TRUE){
        for(att_e in setdiff(network::list.network.attributes(time_slice), c("bipartite", "directed", "hyper", "loops", "mnext", "multiple", "n"))){
          # setting edge attributes harder given how they're stored
          if(att_e != "na"){
            # el <- as.matrix(time_slice,matrix.type="edgelist")
            #el[,1] <- network::get.vertex.attribute(time_slice, 'vertex.names')[el[,1]]
            #el[,2] <- network::get.vertex.attribute(time_slice, 'vertex.names')[as.numeric(el[,2])]
            #el <- cbind(el, network::get.edge.attribute(time_slice, att_e))
            #class_att <- class(network::get.edge.attribute(time_slice, att_e))
            #el_red <- as.matrix(red_net, matrix.type="edgelist")
            #el_red[,1] <- network::get.vertex.attribute(red_net, 'vertex.names')[el_red[,1]]
            #el_red[,2] <- network::get.vertex.attribute(red_net, 'vertex.names')[as.numeric(el_red[,2])]
            #el_red <- cbind(el_red, NA)

            adj <- as.matrix(time_slice, matrix.type = "adjacency")
            adj <- adj[order(as.integer(colnames(adj))),order(as.integer(colnames(adj)))]

            net_att <- network::get.network.attribute(time_slice, att_e)

            colnames(net_att) <- get.vertex.attribute(orig_nets[[t]], 'vertex.names')
            rownames(net_att) <- get.vertex.attribute(orig_nets[[t]], 'vertex.names')

            adj_red <- as.matrix(red_net, matrix.type = "adjacency")
            vertices_red <- rownames(adj_red)
            net_red <- net_att[vertices_red, vertices_red]

            # check to make sure that el's are same order
            # if(all(el_red[,1] %in% el[,2])){
            #  el_red <- cbind(el_red[,2], el_red[,1], NA)
            #}

            # aad a third column to el_red which contains the edge attribute from el when edges were kept
            #  m3 <- rbind(el, el_red)
            # el <- m3[!duplicated(m3[,1:2], fromLast = TRUE), , drop = TRUE]
            network::set.network.attribute(red_net,att_e,net_red)

          }
        }
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
        time_indices <- which(keep_mat[i,] == TRUE)
        index <- time_indices[t]
        indices <- which(network::get.vertex.attribute(orig_nets[[index]], 'vertex.names') %in% vertex_ids)
        for(att in network::list.vertex.attributes(orig_nets[[index]])){
          network::set.vertex.attribute(xt[[i]][[t]], att, network::get.vertex.attribute(orig_nets[[index]], att)[indices])
        }
        if(edge_covariates == TRUE){
          for(att_e in setdiff(network::list.network.attributes(orig_nets[[index]]), c("bipartite", "directed", "hyper", "loops", "mnext", "multiple", "n"))){
            if(att_e != "na"){
              #  el <- as.matrix(net[[t]],matrix.type="edgelist")
              # el[,1] <- network::get.vertex.attribute(net[[t]], 'vertex.names')[el[,1]]
              #  el[,2] <- network::get.vertex.attribute(net[[t]], 'vertex.names')[as.numeric(el[,2])]
              #  el <- cbind(el, network::get.edge.attribute(net[[t]], att_e))
              #  att_class <- class(network::get.edge.attribute(net[[t]], att_e))
              #  el_red <- as.matrix(xt[[i]][[t]], matrix.type="edgelist")
              #  el_red[,1] <- network::get.vertex.attribute(xt[[i]][[t]], 'vertex.names')[el_red[,1]]
              #  el_red[,2] <- network::get.vertex.attribute(xt[[i]][[t]], 'vertex.names')[as.numeric(el_red[,2])]
              #  el_red <- cbind(el_red, NA)

              adj <- as.matrix(orig_nets[[index]], matrix.type = "adjacency")
              net_att <- network::get.network.attribute(orig_nets[[index]], att_e)
              colnames(net_att) <- colnames(adj)[1:nrow(net_att)]
              rownames(net_att) <- rownames(adj)[1:nrow(net_att)]

              adj_red <- as.matrix(xt[[i]][[t]], matrix.type = "adjacency")
              vertices_red <- stats::na.omit(rownames(adj_red))
              net_red <- net_att[vertices_red, vertices_red]

              # if(all(el_red[,1] %in% el[,2])){
              #    el_red <- cbind(el_red[,2], el_red[,1], NA)
              #  }

              #if(network::network.size(xt[[i]][[t]]) == 1){
              #  el_red <- matrix(nrow = 0, ncol = 3)
              #}

              # m3 <- rbind(el, el_red)
              #el <- m3[!duplicated(m3[,1:2], fromLast = TRUE), , drop = TRUE]

              # if(class(el) == "numeric" & length(el) == 3){
              #    el <- t(as.matrix(el))
              #  }

              # if(att_class == "numeric"){
              #    network::set.edge.attribute(xt[[i]][[t]],att_e,as.numeric(el[,3]))
              #  } else {
              #    network::set.edge.attribute(xt[[i]][[t]],att_e,el[,3])
              #  }
              network::set.network.attribute(xt[[i]][[t]],att_e,net_red)

            }

          }

        }

      }
      if(length(nets) == 1 & t == 1){
        vertex_ids <- network::get.vertex.attribute(nets[[1]], 'vertex.names')
        time_indices <- which(keep_mat[i,] == TRUE)
        index <- time_indices[1]
        indices <- which(network::get.vertex.attribute(orig_nets[[index]], 'vertex.names') %in% vertex_ids) #
        for(att in network::list.vertex.attributes(orig_nets[[index]])){
          if(att != "na"){
            network::set.vertex.attribute(nets[[1]], att, network::get.vertex.attribute(orig_nets[[index]], att)[indices])
          }
        }
        if(edge_covariates == TRUE){
          for(att_e in setdiff(network::list.network.attributes(time_slice), c("bipartite", "directed", "hyper", "loops", "mnext", "multiple", "n"))){
            if(att_e != "na"){
              #el <- as.matrix(net[[t]],matrix.type="edgelist")
              #el[,1] <- network::get.vertex.attribute(net[[t]], 'vertex.names')[el[,1]]
              #el[,2] <- network::get.vertex.attribute(net[[t]], 'vertex.names')[as.numeric(el[,2])]
              #el <- cbind(el, network::get.edge.attribute(net[[t]], att_e))
              #
              #attr_class <- class(network::get.edge.attribute(net[[t]], att_e))
              #
              #el_red <- as.matrix(nets[[1]], matrix.type="edgelist")
              #el_red[,1] <- network::get.vertex.attribute(nets[[1]], 'vertex.names')[el_red[,1]]
              #el_red[,2] <- network::get.vertex.attribute(nets[[1]], 'vertex.names')[as.numeric(el_red[,2])]
              #el_red <- cbind(el_red, NA)


              adj <- as.matrix(orig_nets[[index]], matrix.type = "adjacency")
              net_att <- network::get.network.attribute(orig_nets[[index]], att_e)
              colnames(net_att) <- colnames(adj)[1:nrow(net_att)]
              rownames(net_att) <- rownames(adj)[1:nrow(net_att)]

              adj_red <- as.matrix(nets[[1]], matrix.type = "adjacency")
              vertices_red <- stats::na.omit(rownames(adj_red))
              net_red <- net_att[vertices_red, vertices_red]

              #if(all(el_red[,1] %in% el[,2])){
              #  el_red <- cbind(el_red[,2], el_red[,1], NA)
              #}
              #
              #if(network::network.size(xt[[i]][[t]]) == 1){
              #  el_red <- matrix(nrow = 0, ncol = 3)
              #}
              #
              #m3 <- rbind(el, el_red)
              #el <- m3[duplicated(m3[,1:2], fromLast = TRUE), , drop = TRUE]
              #
              #if(class(el) == "numeric" & length(el) == 3){
              #  el <- t(as.matrix(el))
              #}
              #
              #if(attr_class == "numeric"){
              #  network::set.edge.attribute(nets[[1]],att_e,as.numeric(el[,3]))
              #} else {
              #network::set.edge.attribute(nets[[1]],att_e,el[,3])
              #}
              network::set.network.attribute(nets[[1]],att_e,net_red)


            }

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
  if(length(dropped) > 0){
    cat("Longitudinally observed ego-networks removed because no time periods achieve the minimum size necessary to be included and to allow for an identifiable model.  If necessary, adjust min_size parameter.", "\n")
    print(paste(dropped, "removed from ego-TERGM analysis due to no identifiable networks"))
  }

  remaining_vertices <- vertices[-c(null_networks)]

  xt <- xt[!is.na(xt)]


  rm(ego, keep, time_slice, red_net, i, t, red_net_list, nets, vertex_ids, indices)


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

        fit <- btergm_local(form, R=R, parallel=parallel, ncpus=ncpus, verbose = FALSE)
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
      cat("Finished kmeans initialization.", "\n")
    }
    return(list(theta=theta, group.theta=group.theta, lambda=lambda, roles=roles, fit.list=fit.list))
  }


  init<-init.egoergm(form = form, roles = roles, p.ego.terms = NULL, R = R, parallel = parallel, ncpus = ncpus, seed = seed)

  ########################################################################
  ### Calculating and save change statistics
  ########################################################################

  ergmformula <- paste("~", paste(form,collapse="+"),sep="") # Establish function ergm formula that includes the ego.terms object
  obs.S<-list() #obs.S list that will be useful for the below for loop
  cat("Calculating all network statistics.", "\n")
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

  cat("Network statistics calculated.", "\n")

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

  cat("EM algorithm starting.", "\n")

  out <- fit.mix.egoergm(form=form,init=init,obs.S,roles) # run the EM algorithm
  lambda<-out$lambda
  group.theta<-out$theta
  EE.BIC<-out$EE.BIC
  TS.BIC <- out$TS.EE.BIC
  z<-apply(lambda, 1, which.max)
  roles_out <- data.frame(Id = remaining_vertices,
                          Role = z)

  cat("EM algorithm completed.", "\n")
  cat("Done.", "\n")
  cat("Completed Time:", format(Sys.time(), "%a %b %d %X %Y"), "\n")
  return(list(model.fit = "egoTERGM",
              lambda = lambda,
              group.theta = group.theta,
              EE.BIC = EE.BIC,
              TS.BIC = TS.BIC,
              role_assignments = roles_out,
              reduced_networks = reduced_networks,
              form = form,
              ego_nets = x,
              ego_nets_used = keep_mat))

}
