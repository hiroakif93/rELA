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
#' @export
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
        cat(sprintf('Done ; elapsed time %.2f sec\n\n', proc.time()[3]-s))	
        ## ============================================== ##

    }else{
        
        ## ============================================== ##
        ## -- with explicit variables
        fittingMat <- matrix(0, ncol=ncol(data)+2+(ncol(env)-1), nrow=ncol(data),
        		     dimnames=list(colnames(data), c('h', paste('g',1:(1+(ncol(env)-1)), sep='.'), 
								paste('J',colnames(data),sep='.')) ))
        					
        cat('Start parameter fitting\n')	
        s <- proc.time()[3]			
        for(i in 1:rep){
            
            cat(sprintf('Fitting %s...', i))	
            SAres <- SA(ocData = data, envData=as.matrix(env), maxInt = max.itr)
            fittingMat <- fittingMat + SAres
        }
        fit <- fittingMat/rep
        cat(sprintf('Done ; elapsed time %.2f sec\n\n', proc.time()[3]-s))	
        ## ============================================== ##
        
    }
    
    return(fit)
}

#'Parallel stochastic Approximation
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach foreach
#' @export
runSAparallel <- function(data=NULL, env=NULL, 
                   	  rep=16, max.itr=10000, thread=1){
	
	require(doParallel)
	for.parallel(thread)
	
    if(is.null(env)){
    	
    	## ============================================== ##
        ## -- without explicit variables
        
        cat('Start parameter fitting\n')	
        s <- proc.time()[3]			
        fittingMat <- foreach(i = 1:rep) %dopar% {
            
	    cat( sprintf('Trial %s...',i) )
            SAres <- SA_simple(ocData = data, maxInt = max.itr)
            return(SAres)
        }
        
        cat(sprintf('\nDone ; elapsed time %.2f sec\n\n', proc.time()[3]-s))	
        ## ============================================== ##
		fittingRes <- matrix(0, ncol=ncol(data)+1, nrow=ncol(data),
        		     dimnames=list(colnames(data), c('h', paste('J',colnames(data),sep='.')) ))
        					
    	for(i in 1:rep)fittingRes <- fittingRes + fittingMat[[i]]
	    
    }else{
        
        ## ============================================== ##
        ## -- with explicit variables
				
        cat('Start parameter fitting\n')	
        s <- proc.time()[3]			
        fittingMat <- foreach(i = 1:rep) %dopar% {
		
            cat( sprintf('Trial %s...',i) )
            SAres <- SA(ocData = data, envData=as.matrix(env), maxInt = max.itr)
            return(SAres)
        }
      
        cat(sprintf('\nDone ; elapsed time %.2f sec\n\n', proc.time()[3]-s))	
        ## ============================================== ##
	if( ncol(as.matrix(env))==1 ){ 
	    gname <- paste('g', colnames(env), sep='.') 
	    ncols <- ncol(data)+2
	}else{
	    gname <- paste('g',1:(1+(ncol(env)-1)), sep='.')	
	    ncols <- ncol(data)+2+(ncol(env)-1)
	}
	    
        fittingRes <- matrix(0, ncol= ncols, nrow=ncol(data),
        		     dimnames=list(colnames(data), c('h', gname, 
					   paste('J',colnames(data),sep='.')) ))
        					
    	for(i in 1:rep)fittingRes <- fittingRes + fittingMat[[i]]
    }
       					
    return(fittingRes/rep)
}

