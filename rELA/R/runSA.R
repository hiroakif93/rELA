#'Stochastic Approximation
#'@description Estimating species interaction, environemnt preference.
#'
#'@param data : data is a matrix. Row is sample, column is species.
#'@param env : env is a matrix. You should standardize it.
#'@param rep ; rep is a integer.
#'@param max.itr : max.itr is a integer. 
#'
#' @useDynLib rELA, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' 
runSA <- function(data=NULL, env=NULL,
                  rep=16, max.itr=10000){
    
    if(is.null(env)){
    	
    	## ============================================== ##
        ## -- without explicit variables
        fittingMat <- matrix(0, ncol=ncol(data)+1, nrow=ncol(data),
        					dimnames=list(colnames(data), c('h', paste('J',colnames(data),sep='.')) ))
        					
        cat('Start parameter fitting\n')	
        s <- proc.time()[3]			
        for(i in 1:rep){
            
            cat(sprintf('Fitting %s...', i))	
            SAres <- SA_simple(ocData = data, maxInt = max.itr)
            fittingMat <- fittingMat + SAres
        }
        fit <- fittingMat/rep
        cat(sprintf('Done ; elapsed time %.2f sec', proc.time()[3]-s))	
        ## ============================================== ##

    }else{
        
        ## ============================================== ##
        ## -- with explicit variables
        fittingMat <- matrix(0, ncol=ncol(data)+2, nrow=ncol(data),
        					dimnames=list(colnames(data), c('h', 'g', paste('J',colnames(data),sep='.')) ))
        					
        cat('Start parameter fitting\n')	
        s <- proc.time()[3]			
        for(i in 1:rep){
            
            cat(sprintf('Fitting %s...', i))	
            SAres <- SA(ocData = data, envData=as.matrix(env), maxInt = max.itr)
            fittingMat <- fittingMat + SAres
        }
        fit <- fittingMat/rep
        cat(sprintf('Done ; elapsed time %.2f sec', proc.time()[3]-s))	
        ## ============================================== ##
        
    }
    
    return(fit)
}
