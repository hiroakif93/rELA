#'Energy landscape analysis
#'@description Estimating stable state and tipping point
#'
#'@param alpha : alpha is a vector which explicit/implicit preference. Use output of runSA function.
#'@param J : J is a matrix which shows species preference. Use output of runSA function.
#'@param SS.itr : SS.itr is a integer of iterator of stable state estimate.
#'@param FindingTip.itr : SS.itr is a integer of iterator of tipping points estimate.
#'
#' @useDynLib rELA, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @export
ELA <- function(alpha=alpha, J=jj, 
                SS.itr=20000,
                FindingTip.itr=10000){ 
    start <- proc.time()[3]
    ## ================================ ##
    ## -- Stable state
    minsets=SSestimate(alpha, J, itr=SS.itr)
    minsets=unique(round(minsets, digits=8))
    stablestate=t(apply(minsets, 1, 
                        function(x){ ssid=strtoi( paste0(x[-ncol(minsets)], collapse=''), base=2)
                        c(ssid, x[ncol(minsets)]) }))
    
    ## ================================ ##
    ## -- Tipping point
    comb=expand.grid(1:nrow(minsets), 1:nrow(minsets))
    tippointSummary <- tippoint <- c()
    for(k in 1:nrow(comb)){#k=3
        if(comb[k,1] < comb[k,2] ){
            ss1=minsets[comb[k,1], -ncol(minsets)]
            ss2=minsets[comb[k,2], -ncol(minsets)]
                                                              
            tippoint.tmp=FindingTippingpoint_cpp(s1=ss1, s2=ss2,
                                          alpha=alpha, jj=J, 
                                          tmax=FindingTip.itr)
            tipstate <- t(apply(tippoint.tmp, 1, 
                                function(x){ tipid=strtoi( paste0(x[-ncol(minsets)], collapse=''), base=2)
                                c(tipid, x[ncol(minsets)]) }))
            tippointSummary <- rbind(tippointSummary, tipstate)
            tippoint <- rbind(tippoint, tippoint.tmp)
        }else{
            tippointSummary <- rbind(tippointSummary, c(Inf, Inf))
        }
    }
    
    rbind(data.frame('State'='stable', stablestate), data.frame('State'='tip', tippointSummary))
    ## ================================ ##
    end <- proc.time()[3]
    cat(sprintf('Elapsed time %.2f sec\n', end-start))
}


#'Calculation stability indices.
#'@description Claculating stable state energy, difference between sample's energy and stable state energy and stable state entropy.
#'
#'@param data : data is a matrix.
#'@param alpha : alpha is a vector which explicit/implicit preference. Use output of runSA function.
#'@param J : J is a matrix which shows species preference. Use output of runSA function.
#'@param 
#'
#' @useDynLib rELA, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @export
calcStability <- function(data=NULL, alpha, J,
						  seitr=1000, convTime=10000){
						  	
    start <- proc.time()[3]
    ## ================================ ##
    sampleSS <- t(apply(data, 1, SteepestDescent_cpp, alpha=alpha, beta=J))
    ssid <- apply(sampleSS[,-ncol(sampleSS)], 1, 
                  function(x){ paste0(x, collapse='') })
    ssenergy <- sampleSS[, ncol(sampleSS)]
    energy <- apply(ocdata, 1, cEnergy, alpha= alpha, beta=J)
    
    energy.gap <- energy-ssenergy
    ssent <- SSentropy(uoc=ocdata, ss=unique(sampleSS[,-ncol(sampleSS)]),
          				alpha= alpha, beta=J, 
          				seitr=seitr, convTime=convTime)
          	  
    return( data.frame(StableState.id=ssid, SSenergy=ssenergy,
                       energy.gap=energy.gap, SSentropy= ssent[,1]))
                       
    ## ================================ ##
    end <- proc.time()[3]
    cat(sprintf('Stability calculation done.\nElapsed time %.2f sec\n\n', end-start))                   
}