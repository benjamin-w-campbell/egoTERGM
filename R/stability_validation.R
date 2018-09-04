#' Stability validation of ego-Temporal Exponential Random Graph Model (ego-TERGM) fit.
#'
#' This function examines the stability of ego-TERGM results to the pooling window.  It takes some proportion of the window's network and compares the result of a model fit on these time periods to the original fit.
#' @param ego_tergm_fit The output of a previously fit ego-TERGM fit using the ego_tergm function.  This is the model that stability validation will be performed on.
#' @param splitting_probability A value from 0 to 1 that determines the probability that any given network is assigned to the comparison group.
#' @param seed The seed set to replicate analysis for pseudorandom number generator.
#' @return Returns a matrix of cross-tabulation results that allows for comparisons of common cluster assignments.  Note that labels may switch.
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
#'                           parallel = "no",
#'                           ncpus = 1,
#'                           steps = 50,
#'                           tol = 1e-06)
#' }
#' @export
