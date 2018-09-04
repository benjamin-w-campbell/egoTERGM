#' Custom interpret function for ego-TERGM.
#'
#' This function assists the user in interpreting output from the ego_tergm function.
#' @param ego_tergm_fit The output from a fitted "ego_tergm".
#' @param custom_var_names A vector of character terms in the same order as the form object fed to ego_tergm of clearer names for these variables.
#' @keywords summary interpretation
#' @references{
#' Campbell, Benjamin W. (2018):
#'  Inferring Latent Roles in Longitudinal Networks.
#'  \emph{Political Analysis} 26(3): 292-311.  \url{https://doi.org/10.1017/pan.2018.20}
#' }
#' @examples
#' interpret_ego_tergm(ego_tergm_fit = ego_tergm_fit)
#' interpret_ego_tergm(ego_tergm_fit = ego_tergm_fit,
#'                     custom_var_names = c("Edges", "Mutual", "Triangle",
#'                                         "In-Degree", "Out-Degree", "Sex Homophily"))
#' @export

interpret_ego_tergm <- function(ego_tergm_fit = NULL, custom_var_names = NULL){
  role_names <- paste0("Role", rep(1:ncol(ego_tergm_fit$lambda)))
  mixing_proportions <- as.matrix(lapply(c(1:3), function(x) mean(ego_tergm_fit$lambda[,x])), nrow = ncol(ego_tergm_fit$lambda), ncol = 1)

  if(is.null(custom_var_names)){
    centroids <- cbind(mixing_proportions, ego_tergm_fit$group.theta)
    colnames(centroids) <- c("Mixing Proportions", ego_tergm_fit$form)
    rownames(centroids) <- role_names
  } else {
    centroids <- cbind(mixing_proportions, ego_tergm_fit$group.theta)
    colnames(centroids) <- c("Mixing Proportions", custom_var_names)
    rownames(centroids) <- role_names
  }
  cat("EM role-based centroids.  Note: Do not interpret as unbiased ERGM parameters, see Salter-Townshend and Murphy (2015) and Campbell (2018) for details.\n")
  print(centroids)
}
