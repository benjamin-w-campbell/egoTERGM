#' Simulation of overlapping longitudinal ego-networks according to a mixture of data generating processes and fitting of ego-TERGM to that network.
#'
#' This function simulated longitudinal ego-networks that are dependent and overlap, and estimates an ego-TERGM. Useful for monte carlos and proofs of concept.
#' @param form A vector of ERGM terms that should be used to generate the networks.  A homophily term will be added to this formula to estimate the ego-TERGM.
#' @param params A "roles" by "form" matrix of ERGM simulation parameters with groups defined on the row and model terms defined on the column.
#' @param roles An integer for the number of distinct mixture groups that should be simulated.
#' @param N_per_role An integer for the number of different longitudinally observed ego-networks that should be simulated per role in roles.
#' @param t_steps An integer for the number of time steps that each ego-network should be observed across.
#' @param egonet_size An integer for the size of each ego-network simulated.
#' @param seed The seed set to replicate analysis for pseudorandom number generator.
#' @param prop_fixed A proportion (0,1 bound) of ties in overlapping ego-networks that should be forced "on".  For overlapping ego-networks this proportion of ties is assumed to follow a homophilous tie-forming proceedure.  A homophily
#' @param N_overlap An integer for the number of ego-networks for a role that should overlap with each data-generating process (defined by roles-1).
#' @param R The number of bootstrap replications that should be used for the estimation of a bootstrapped MPLE estimated TERGM for model initialization.  Defaults to 10.
#' @param parallel How the BTERGM should be computed, either with parallel processing ("multicore", "snow") or no parallel processing ("no").  Defaults to no.
#' @param ncpus The number of CPUs that should should be used for estimation, defaults to 1.
#' @param steps The number of default EM steps that should be taken, defaults to 50.
#' @param tol The difference in parameter estimates between EM iterations to determine if the algorithm has converged.  Defaults to 1e-6.
#' @param seed The seed for simulation
#' @return A list of simulated ego-networks and the output of the ego_tergm function fit to this.
#' @keywords simulation networks
#'# @references Campbell, Benjamin W. 2017. Inferring Latent Roles in Longitudinal Networks using the Ego-TERGM. Working Paper.
#' @examples
#' \dontrun{
#' net <- sim_overlapping_egonets(form = c("edges", "gwesp(0.8,fixed=TRUE)", "gwdegree(decay=0.8,fixed=TRUE)"),
#'                   params = rbind(c(-3,1,0), c(-1,-2,-1), c(-2,0,2)),
#'                   roles = 3,
#'                   N_per_role = 10,
#'                   t_steps = 3,
#'                   egonet_size = 20,
#'                   R = 30,
#'                   parallel = "multicore",
#'                   ncpus = 3,
#'                   steps = 50,
#'                   tol = 1e-6)
#' }
#' @export
sim_overlap_egonets <- function(form = NULL, params = NULL, roles = NULL, N_per_role = NULL, t_steps = NULL, egonet_size = NULL,
                        R = 10, parallel = "no", ncpus = 1, steps = 50, tol = 1e-6, seed = 12345){
  set.seed(seed)


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
      sim.net[[i]]<- simulate.formula(as.formula(paste("network(NN,directed=FALSE)", ergmformula)), coef = c(params[(1:roles)[sim.K[i]],]))
    }
    sim.x[[k]] <- sim.net
  }
  sim.K <- sim.x





}
