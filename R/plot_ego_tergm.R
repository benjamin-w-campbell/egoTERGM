#' Function to plot role assignments.
#'
#' This function assists the user in interpreting output from the ego_tergm function.
#' @param ego_tergm_fit The output from a fitted "ego_tergm".
#' @param plot_indices A vector of indices reflecting the time steps that plots should be returned for.  In networks observed over many time intervals it may be unrealistic to plot every network.
#' @param node.size The size of the nodes fed to ggnet2.
#' @param edge.size The size of the edges fed to ggnet2.
#' @param edge.color The color of edges fed to ggnet2.
#' @keywords summary interpretation plot assignments
#' @return A list of network plots with nodes colored by role assignment.
#' @references{
#' Campbell, Benjamin W. (2018):
#'  Inferring Latent Roles in Longitudinal Networks.
#'  \emph{Political Analysis} 26(3): 292-311.  \url{https://doi.org/10.1017/pan.2018.20}
#' }
#' @examples
#' plots <- plot_ego_tergm(ego_tergm_fit = ego_tergm_fit)
#' plots[[1]]
#' plots[[2]]
#' @export

plot_ego_tergm <- function(ego_tergm_fit = NULL, plot_indices = NULL, node.size = 6, edge.size = 1, edge.color = "grey"){

  plot_function <- function(nw = NULL){
    if(is.network(nw)){
      p <- GGally::ggnet2(nw, color = "role", label = TRUE, node.size = node.size, edge.size = edge.size, edge.color = edge.color)
    } else {
      p <- NA
    }
    return(p)
  }

  append_roles <- function(nw = NULL){
    nw <- nw[[1]]
    if(is.network(nw)){
      vertices <- data.frame(Id = as.character(get.vertex.attribute(nw, "vertex.names")))
      role_df <- merge(vertices, ego_tergm_fit$role_assignments, by = "Id", all.x = TRUE, all.y = FALSE)
      role_df <- role_df[order(match(role_df[,1],vertices[,1])),]
      net <- set.vertex.attribute(x = nw, attrname = "role", value = role_df$Role)
    } else {
      net <- NA
    }
    return(net)
  }

  if(is.null(plot_indices)){
    plot_indices <- rep(1:length(ego_tergm_fit$reduced_networks))
  }

  # append the role vector to each network for colors as "role"
  net_list <- list()
  for(i in 1:length(ego_tergm_fit$reduced_networks)){
    temp <- append_roles(ego_tergm_fit$reduced_networks[i])
    net_list[[i]] <- temp
  }

  # create plot_list
  plot_list <- lapply(net_list, function(x) plot_function(x))

   return(plot_list)
}
