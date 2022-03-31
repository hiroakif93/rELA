#'Registering prallel
#' @export
#' 
for.parallel <- function (cores = 1, clusterType='FORK') 
{
    require(doParallel)
    cluster = makeCluster(cores, clusterType)
    return(registerDoParallel(cluster))
}