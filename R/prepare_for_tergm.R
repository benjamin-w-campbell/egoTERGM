#' Prepares ego-TERGM output for xergm's btergm function.
#'
#' This takes the output of an ego-TERGM call and prepares it for use by the xergm btergm function.  Note: This routine assumes temporal independence within ego-networks and independence across ego-networks.
#' @param ego_tergm_fit The output of an ego-TERGM call.
#' @return A list of length G containing pooled cluster assignments.  First-level elements of this list may be fed to a btergm call.
#' @references{
#' Campbell, Benjamin W. (2018):
#'  Inferring Latent Roles in Longitudinal Networks.
#'  \emph{Political Analysis} 26(3): 292-311.  \url{https://doi.org/10.1017/pan.2018.20}
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
#'                           forking = FALSE,
#'                           ncpus = 1,
#'                           steps = 50,
#'                           tol = 1e-06)
#'
#' net_list <- prepare_for_tergm(ego_tergm_fit)
#'
#' role1_btergm <- btergm(net_list[[1]] ~ edges + mutual + triangle + nodeicov('idegsqrt') +
#'                                        nodeocov('odegsqrt') + nodematch('sex'),
#'                        R = 500)
#'
#' }
#' @export

prepare_for_tergm <- function(ego_tergm_fit = NULL){
  cat("Note: The output of this function cannot be used in a TERGM containing memory terms and assumes temporal independence within ego-networks and independence across ego-networks.")
  tergm_list <- list()
  G_detect <- length(unique(ego_tergm_fit$role_assignments$Role))
  roles <- unique(ego_tergm_fit$role_assignments$Role)

  flattenlist <- function(x){
    morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
    out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
    if(sum(morelists)){
      Recall(out)
    } else {
      return(out)
    }
  }


  for(i in 1:G_detect){
    indices <- which(ego_tergm_fit$role_assignments$Role == roles[i])
    ego_nets_roles <- ego_tergm_fit$ego_nets[indices]
    pooled_networks <- flattenlist(ego_nets_roles)
    tergm_list[[i]] <- pooled_networks
  }

  return(tergm_list)
}
