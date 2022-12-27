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
	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Stable state estimatin	
	minsets = SSestimate(alpha, J, itr = SS.itr)
	minsets = unique(minsets)

	uniqueSS <- data.frame(unique(minsets[,-ncol(minsets)]))
	uniqueSS$energy <- NA

	for(i in 1:nrow(uniqueSS)){

	    state <- uniqueSS[i,-ncol(uniqueSS)]
	    detect <- apply(minsets[,-ncol(minsets)], 1, function(x){ all(x==state) })

	    uniqueSS[i, ncol(uniqueSS)] <- mean(minsets[detect, ncol(minsets)])
	}
	sampleSS <- uniqueSS[order(uniqueSS$energy),]

	minsets <- as.matrix(sampleSS[order(sampleSS[, ncol(sampleSS)]), ])
	colnames(minsets) <- c(speciesName, "energy")

	ssid <- sprintf("SS_%s", formatC(1:nrow(minsets), width = nchar(ncol(minsets)), 
					 flag = "0"))
	rownames(minsets) <- ssid

	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Tipping point
	comb=expand.grid(1:nrow(minsets), 1:nrow(minsets))
	comb <- comb[comb[,1]!=comb[,2],2:1]

	tpmat <- matrix(Inf, ncol=nrow(comb), nrow=nrow(minsets))
	dimnames(tpmat) <- list(rownames(minsets), 1:ncol(tpmat))

	tpnodeID <- sprintf("TPnode_%s", formatC(1:nrow(comb), width = nchar(nrow(comb)), 
						 flag = "0"))

	FindTip <- c()
	for(k in 1:nrow(comb)){

	    minsetsub <- minsets[as.integer(comb[k, ]), ]
	    ss1 = minsetsub[1, -ncol(minsets)]
	    ss2 = minsetsub[2, -ncol(minsets)]

	    tippoint.tmp = FindingTippingpoint_cpp(s1 = ss1,  s2 = ss2, 
						   alpha = alpha, jj = J,
						   tmax = FindingTip.itr)

	    colnames(tippoint.tmp) <- c(speciesName, "energy")

	    ss <- data.frame(SS1=rownames(minsetsub)[1], SS1.energy=minsetsub[1,ncol(minsetsub)],
			     SS2=rownames(minsetsub)[2], SS2.energy=minsetsub[2,ncol(minsetsub)])
	    tp <- data.frame(TP=tpnodeID[k], tippoint.tmp)
	    FindTip <- rbind(FindTip, cbind(ss, tp))

	}

	## ================================================ ##
	## -- Summarize
	tpState <- (FindTip[,-c(1:5)]);
	tpStateUni <- unique(tpState)
	tpid <- sprintf("TP_%s", formatC(1:nrow(tpStateUni), width = nchar(nrow(tpStateUni)), 
					 flag = "0"))
	rownames(tpStateUni) <- tpid

	stateInfo <- rbind(data.frame(stateID=rownames(minsets), state=rep('Stable state', nrow(minsets)), minsets), 
			   data.frame(stateID=rownames(tpStateUni), state=rep('Tipping point', nrow(tpStateUni)), tpStateUni))


	# ||||||||||||||||||||||||||||||||||||| ##
	Stable <- FindTip[,c(1:4)]

	tmp <- rownames(tpStateUni); names(tmp) <-apply(tpStateUni[,-ncol(tpStateUni)], 1, paste, collapse='')
	ordered <- apply(tpState[,-ncol(tpState)], 1, paste, collapse='')
	tpidVec <- tmp[ordered]

	network <- data.frame(Stable, TP=tpidVec, Tp.energy=FindTip[,ncol(FindTip)])
	elasummary <- list(stateInfo = stateInfo, network =network)

	# ||||||||||||||||||||||||||||||||||||| ##
	end <- proc.time()[3]
	cat(sprintf("Elapsed time %.2f sec\n", end - start))

	return(elasummary)
}

# ELA <- function(alpha=alpha, J=jj, 
#                 SS.itr=20000,
#                 FindingTip.itr=10000){ 
#     start <- proc.time()[3]
#     
#     speciesName <- rownames(J)
#     
#     ## ||||||||||||||||||||||||||||||||||||| ##
#     ## -- Stable state estimatin	
#     minsets = SSestimate(alpha, J, itr = SS.itr)
#     minsets = unique(minsets)
#     
#     uniqueSS <- data.frame(unique(minsets[,-ncol(minsets)]))
#     uniqueSS$energy <- NA
#     
#     for(i in 1:nrow(uniqueSS)){
#         
#         state <- uniqueSS[i,-ncol(uniqueSS)]
#         detect <- apply(minsets[,-ncol(minsets)], 1, function(x){ all(x==state) })
#         
#         uniqueSS[i, ncol(uniqueSS)] <- mean(minsets[detect, ncol(minsets)])
#     }
#     sampleSS <- uniqueSS[order(uniqueSS$energy),]
#     
#     minsets <- as.matrix(sampleSS[order(sampleSS[, ncol(sampleSS)]), ])
#     colnames(minsets) <- c(speciesName, "energy")
#     
#     ssid <- sprintf("SS_%s", formatC(1:nrow(minsets), width = nchar(ncol(minsets)), 
#                                      flag = "0"))
#     rownames(minsets) <- ssid
#     
#     ## ||||||||||||||||||||||||||||||||||||| ##
#     ## -- Tipping point estimation
#     
#     comb = expand.grid(1:nrow(minsets), 1:nrow(minsets))
#     comb <- comb[comb[, 1] <comb[, 2], 2:1]
#     
#     tpnodeID <- sprintf("TPnode_%s", formatC(1:nrow(comb), width = nchar(nrow(comb)), 
#                                              flag = "0"))
#     
#     
#     FindTipRes <- TPestimate(as.matrix(comb), minset = minsets[,-ncol(minsets)],
#                              alpha, J, 100000)
#     
#     FindTip <- data.frame(matrix(NA, ncol=4, nrow=nrow(FindTipRes)) ,TP=tpnodeID, FindTipRes)
#     colnames(FindTip) <- c('SS1', 'SS1.energy', 'SS2', 'SS2.energy','TP',speciesName, "energy")
#     
#     for (k in 1:nrow(comb)) {
#         minsetsub <- minsets[as.integer(comb[k, ]), ]
#         ss1 = minsetsub[1, -ncol(minsets)]
#         ss2 = minsetsub[2, -ncol(minsets)]
#         
#         
#         ss <- data.frame(SS1=rownames(minsetsub)[1], SS1.energy=minsetsub[1,ncol(minsetsub)],
#                          SS2=rownames(minsetsub)[2], SS2.energy=minsetsub[2,ncol(minsetsub)])
#         FindTip[k, 1:4] <- ss
#     }
#     ## ||||||||||||||||||||||||||||||||||||| ##
#     ## -- Summarize
#     tpState <- (FindTip[,-c(1:5)]);
#     tpStateUni <- unique(tpState)
#     tpid <- sprintf("TP_%s", formatC(1:nrow(tpStateUni), width = nchar(nrow(tpStateUni)), 
#                                      flag = "0"))
#     rownames(tpStateUni) <- tpid
#     
#     stateInfo <- rbind(data.frame(stateID=rownames(minsets), state=rep('Stable state', nrow(minsets)), minsets), 
#                        data.frame(stateID=rownames(tpStateUni), state=rep('Tipping point', nrow(tpStateUni)), tpStateUni))
#     
#     
#     # ||||||||||||||||||||||||||||||||||||| ##
#     Stable <- FindTip[,c(1:4)]
#     
#     tmp <- rownames(tpStateUni); names(tmp) <-apply(tpStateUni[,-ncol(tpStateUni)], 1, paste, collapse='')
#     ordered <- apply(tpState[,-ncol(tpState)], 1, paste, collapse='')
#     tpidVec <- tmp[ordered]
#     
#     network <- data.frame(Stable, TP=tpidVec, Tp.energy=FindTip[,ncol(FindTip)])
#     elasummary <- list(stateInfo = stateInfo, network =network)
#     
#     # ||||||||||||||||||||||||||||||||||||| ##
#     
#     end <- proc.time()[3]
#     cat(sprintf("Elapsed time %.2f sec\n", end - start))
#     
#     return(elasummary)
# }



#'Energy landscape analysis in parallel mode
#'@description Estimating stable state and tipping point
#'
#'@param alpha : alpha is a vector which explicit/implicit preference. Use output of runSA function.
#'@param J : J is a matrix which shows species preference. Use output of runSA function.
#'@param SS.itr : SS.itr is a integer of iterator of stable state estimate.
#'@param FindingTip.itr : SS.itr is a integer of iterator of tipping points estimate.
#'
#' @useDynLib rELA, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach foreach
#' @export
ELAparallel <- function(alpha=alpha, J=jj, 
                SS.itr=20000, FindingTip.itr=10000,
                thread=1){ 
    
    ## ||||||||||||||||||||||||||||||||||||| ##           	
               	
    start <- proc.time()[3]

	speciesName <- rownames(J)
	
	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Stable state estimatin	
	minsets = SSestimate(alpha, J, itr = SS.itr)
	minsets = unique(minsets)
	
	uniqueSS <- data.frame(unique(minsets[,-ncol(minsets)]))
	uniqueSS$energy <- NA
	
	for(i in 1:nrow(uniqueSS)){
	    
	    state <- uniqueSS[i,-ncol(uniqueSS)]
	    detect <- apply(minsets[,-ncol(minsets)], 1, function(x){ all(x==state) })
	    
	    uniqueSS[i, ncol(uniqueSS)] <- mean(minsets[detect, ncol(minsets)])
	}
	sampleSS <- uniqueSS[order(uniqueSS$energy),]
	
	minsets <- as.matrix(sampleSS[order(sampleSS[, ncol(sampleSS)]), ])
	colnames(minsets) <- c(speciesName, "energy")
	
	ssid <- sprintf("SS_%s", formatC(1:nrow(minsets), width = nchar(ncol(minsets)), 
	                                 flag = "0"))
	rownames(minsets) <- ssid
	
	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Tipping point estimation
	
	comb = expand.grid(1:nrow(minsets), 1:nrow(minsets))
	comb <- comb[comb[, 1] <comb[, 2], 2:1]
	
	tpnodeID <- sprintf("TPnode_%s", formatC(1:nrow(comb), width = nchar(nrow(comb)), 
	                             flag = "0"))
	
	for.parallel(thread)
	FindTip <- foreach (k = 1:nrow(comb), .combine=rbind)%dopar%{
	    
	    minsetsub <- minsets[as.integer(comb[k, ]), ]
	    ss1 = minsetsub[1, -ncol(minsets)]
	    ss2 = minsetsub[2, -ncol(minsets)]
	    
	    tippoint.tmp = FindingTippingpoint_cpp(s1 = ss1,  s2 = ss2, 
	                                           alpha = alpha, jj = J,
	                                           tmax = FindingTip.itr)
	
	    colnames(tippoint.tmp) <- c(speciesName, "energy")
	    
	    ss <- data.frame(SS1=rownames(minsetsub)[1], SS1.energy=minsetsub[1,ncol(minsetsub)],
	               SS2=rownames(minsetsub)[2], SS2.energy=minsetsub[2,ncol(minsetsub)])
	    tp <- data.frame(TP=tpnodeID[k], tippoint.tmp)
	    return( cbind(ss, tp))
	
	}
	
	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Summarize
	tpState <- (FindTip[,-c(1:5)]);
	tpStateUni <- unique(tpState)
	tpid <- sprintf("TP_%s", formatC(1:nrow(tpStateUni), width = nchar(nrow(tpStateUni)), 
	                                 flag = "0"))
	rownames(tpStateUni) <- tpid
		
	stateInfo <- rbind(data.frame(stateID=rownames(minsets), state=rep('Stable state', nrow(minsets)), minsets), 
	                   data.frame(stateID=rownames(tpStateUni), state=rep('Tipping point', nrow(tpStateUni)), tpStateUni))
	
	
	# ||||||||||||||||||||||||||||||||||||| ##
	Stable <- FindTip[,c(1:4)]
	
	tmp <- rownames(tpStateUni); names(tmp) <-apply(tpStateUni[,-ncol(tpStateUni)], 1, paste, collapse='')
	ordered <- apply(tpState[,-ncol(tpState)], 1, paste, collapse='')
	tpidVec <- tmp[ordered]
	
	network <- data.frame(Stable, TP=tpidVec, Tp.energy=FindTip[,ncol(FindTip)])
	elasummary <- list(stateInfo = stateInfo, network =network)
	
	# ||||||||||||||||||||||||||||||||||||| ##
	end <- proc.time()[3]
	cat(sprintf("Elapsed time %.2f sec\n", end - start))
	
	return(elasummary)
}



#'ELA Pluning
#'@description Pluning shallow basin
#'
#'@param elasummary : elasummary is a output of ELA function.
#'@param th : th is a threshold which the range is 0 to 1.
#'@export
ELplunning <- function(elasummary, th=0.1){
    
	stateInfo <- elasummary[[1]]
	elaMat <- elasummary[[2]]
	
	## ========================================== ##
	## -- Make incidence matrix
	
	IM <- matrix(Inf, nrow=sum(stateInfo$state=='Stable state'), ncol=nrow(elaMat)+1,
	             dimnames=list(stateInfo$stateID[which(stateInfo$state=='Stable state')],
	                           c('SSenergy', elaMat[,5])))
	IM[,1] <- stateInfo$energy[which(stateInfo$state=='Stable state')]
	for(i in 1:nrow(elaMat)){
	    
	    ss=as.character(elaMat[i,c(1,3)])
	    IM[ss, i+1] <- elaMat[i,ncol(elaMat)]
	    
	}
	
	## ========================================== ##
	## -- Difference between stable state energy
	sub <- IM[,-1]
	eneDiff <- sub-IM[,1]
	
	## ========================================== ##
	## -- Ranking between two stable states
	ssene <- IM[,1]
	pair <-  t(combn(rownames(IM), 2))
	
	pax <- c()
	for(i in 1:nrow(pair)){
	    pairSort <- sort( ssene[pair[i,]])
	    
	    pairEneDiff <- eneDiff[names(pairSort),]
	    tpid <- which(colSums(is.infinite(pairEneDiff))==0)
	    paxtmp <- pairEneDiff[,tpid]
	    
	    pax <- rbind(pax, data.frame(deep=names(paxtmp)[which.max(paxtmp)], 
	                                 shallow=names(paxtmp)[which.min(paxtmp)],  
	                                 min=min(paxtmp), tp=tpid) )
	}
	
	paxmax <- max(pax[,3])
	
	## ========================================== ##
	elaMattmp <- IM
	log <- data.frame(row.names=names(ssene), change=names(ssene))
	t=0
	while( min(pax[,3]) < paxmax*th & nrow(elaMattmp)>2){
	    t=t+1
	    ## |||||||||||||||||||||||||||||||||||||| ##
	    ssenetmp <- elaMattmp[,1]
	    sub <- elaMattmp[,-1]
	    eneDiff <- apply(sub, 2, function(x){x-elaMattmp[,1]})
	    
	    pair <-  t(combn(rownames(elaMattmp), 2))
	    
	    pax <- c()
	    for(i in 1:nrow(pair)){
	        pairSort <- sort( ssenetmp[pair[i,]])
	        
	        pairEneDiff <- eneDiff[names(pairSort),]
	        tpid <- which(colSums(is.infinite(pairEneDiff))==0)
	        paxtmp <- pairEneDiff[,tpid]
	        
	        pax <- rbind(pax, data.frame(deep=names(paxtmp)[which.max(paxtmp)], 
	                                     shallow=names(paxtmp)[which.min(paxtmp)],  
	                                     min=min(paxtmp), tp=tpid) )
	    }
	    
	    paxmax <- max(pax[,3])
	    pp <- pax[which.min(pax[,3]), ]
	    
	    ## |||||||||||||||||||||||||||||||||||||| ##
	    changedState <- log[,ncol(log)]
	    changedState[changedState==as.character(pp['shallow'])] <- as.character(pp['deep'])
	    log <- cbind(log, changedState)
	    
	    elaMattmp <- elaMattmp[-which(rownames(elaMattmp) == as.character(pp['shallow'])),]
	    elaMattmp <- elaMattmp[,-c(as.numeric(pp['tp'])+1)]
	    
	    countInf <- colSums(!is.infinite(elaMattmp))
	    elaMattmp <- elaMattmp[, -which(countInf<2)]
	    
	    ## |||||||||||||||||||||||||||||||||||||| ##
	    
	}
	
	stateid <- unique(c(rownames(elaMattmp), colnames(elaMattmp)[-1]))
	Pluned.stateInfo <- stateInfo[stateid, ]
	
	return(list(network=elaMattmp, stateInfo=Pluned.stateInfo,
	            log=log[,c(1,ncol(log))]))
	## ========================================== ##
    
    
}


#'SteepestDescent
#' @export
SteepestDescent <- function(state, alpha, beta){
	res <- SteepestDescent_cpp(state, alpha, beta)
	return(res)
}

#'cEnergy
#' @export
Energy <- function(state, alpha, beta){
	res <- cEnergy(state, alpha, beta)
	return(res)
}

#'cEnergy
#' @export
SSentropy <- function(state, ss,
          				alpha, beta, 
          				seitr=1000, convTime=10000){
	res <- SSentropy_cpp(uoc= state, ss= ss,
          				alpha= alpha, beta=beta, 
          				seitr=seitr, convTime=convTime)
	return(res)
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
						  plun=0,
						  seitr=1000, convTime=10000){
						  	
    start <- proc.time()[3]
    sampleSS <- t(apply(data, 1, SteepestDescent_cpp, alpha=alpha, beta=J))
   	
   	uniqueSS <- data.frame(unique(sampleSS[,-ncol(sampleSS)]))
	uniqueSS$energy <- NA
	
	for(i in 1:nrow(uniqueSS)){
	    
	    state <- uniqueSS[i,-ncol(uniqueSS)]
	    detect <- apply(sampleSS[,-ncol(sampleSS)], 1, function(x){ all(x==state) })
	    
	    uniqueSS[i, ncol(uniqueSS)] <- mean(sampleSS[detect, ncol(sampleSS)])
	}
	sampleSS <- uniqueSS[order(uniqueSS$energy),]
	
    ## ================================ ##
	if(plun>0){
		
		## |||||||||||||||||||||||||||| ##
		## plunned landscape
	    elasummary <- ELAparallel(alpha=res[,1], J=res[,-1], thread=8)
		redela <- ELplunning(elasummary, th=0.1)
		
		stateInfo <- elasummary[[1]]
		ssInfo <- stateInfo[stateInfo$state=='Stable state',]
		
		## |||||||||||||||||||||||||||| ##
		## -- Convert to deeper landscape
		
		sampleSSpluned <- as.matrix(sampleSS)
		for(i in 1:nrow(ssInfo)){
		    
		    state <- ssInfo[i,-c(1, 2, ncol(ssInfo))]
		    detect <- apply(sampleSS[,-ncol(sampleSS)], 1, function(x){ all(x==state) })
		    
		    if(sum(detect)>1){
		        id <- ssInfo[i,1]
		        
		        if( !any(redela$log[,2]==id)){
		            log <- redela$log
		            change <- log[id,2]
		            
		            sampleSSpluned[detect,] <- matrix( unlist(ssInfo[change, -c(1,2)]), 
		            								   ncol=ncol(sampleSStmp), nrow=sum(detect), 
		            								   byrow=TRUE)
		        }
		        
		    }
		
		}
		sampleSS <- sampleSSpluned
        ## |||||||||||||||||||||||||||| ##
	}
    
    stablestate=t(apply(sampleSS, 1, 
                    function(x){ 
                        sep <- round(length(x)/20)
                        line <- c()
                        for(i in 1:sep){
                            y <- x[-length(x)]
                            bi <- na.omit(y[((1+20*(i-1)):(20*i))])
                            line <- c(line, strtoi( paste0(bi, collapse=''), base=2))
                            
                        }
                        
                        ssid=paste(line, collapse='-')
                        c(ssid, x[length(x)]) }))
                  
    ssenergy <- sampleSS[, ncol(sampleSS)]
    energy <- apply(data, 1, cEnergy, alpha= alpha, beta=J)
    
    energy.gap <- energy-ssenergy
    ssent <- SSentropy_cpp(uoc=data, ss=unique(sampleSS[,-ncol(sampleSS)]),
          				alpha= alpha, beta=J, 
          				seitr=seitr, convTime=convTime)
          	  
    return( data.frame(StableState.id=stablestate[,1], SSenergy=ssenergy,
                       energy.gap=energy.gap, SSentropy= ssent[,1]))
                       
    ## ================================ ##
    end <- proc.time()[3]
    cat(sprintf('Stability calculation done.\nElapsed time %.2f sec\n\n', end-start))                   
}

#'ELAall 
#'@description Run all ELA pipeline
#'@param data : data is a binary matrix.
#'@param alpha : alpha is a vector which explicit/implicit preference. Use output of runSA function.
#'@param J : J is a matrix which shows species preference. Use output of runSA function.
#'@param 
#'
#' @useDynLib rELA, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @export
ELAall <- function( data=NULL, alpha, J,
				    SS.itr=20000, FindingTip.itr=10000,
					plun=0,
					seitr=1000, convTime=10000,
					thread=1){
  start <- proc.time()[3]
    
   elasummary <- ELAparallel(alpha=res[,1], J=res[,-1], 
   			     SS.itr= SS.itr, FindingTip.itr=10000, thread=thread)
   	
    ## ================================ ##    
    sampleSS <- t(apply(data, 1, SteepestDescent_cpp, alpha=alpha, beta=J))
   	
	if(plun>0){
		
		## |||||||||||||||||||||||||||| ##
		## plunned landscape
	    
		redela <- ELplunning(elasummary, th= plun)
		
		stateInfo <- elasummary[[1]]
		ssInfo <- stateInfo[stateInfo$state=='Stable state',]
		
		## |||||||||||||||||||||||||||| ##
		## -- Convert to deeper landscape
		
		sampleSSpluned <- as.matrix(sampleSS)
		ssid <- rep(NA, nrow(sampleSS))
		for(i in 1:nrow(ssInfo)){
		    
		    state <- ssInfo[i,-c(1, 2, ncol(ssInfo))]
		    detect <- apply(sampleSS[,-ncol(sampleSS)], 1, function(x){ all(x==state) })
		    ssid[detect] <- rownames(ssInfo)[i]
        
            if(sum(detect)>1){
		        id <- ssInfo[i,1]
		        
		        if( !any(redela$log[,2]==id)){
		            log <- redela$log
		            change <- log[id,2]
		            
		            sampleSSpluned[detect,] <- matrix( unlist(ssInfo[change, -c(1,2)]), 
		            								   ncol=ncol(sampleSSpluned), nrow=sum(detect), 
		            								   byrow=TRUE)
		        }
		        
		    }
		
		}
		sampleSS <- sampleSSpluned
        ## |||||||||||||||||||||||||||| ##
	}
                      
    ssenergy <- sampleSS[, ncol(sampleSS)]
    energy <- apply(data, 1, cEnergy, alpha= alpha, beta=J)
    
    energy.gap <- energy-ssenergy
    ssent <- SSentropy_cpp(uoc=data, ss=unique(sampleSS[,-ncol(sampleSS)]),
          				alpha= alpha, beta=J, 
          				seitr=seitr, convTime=convTime)
          	  
    stability=data.frame(StableState.id=ssid, SSenergy=ssenergy,
                     energy.gap=energy.gap, SSentropy= ssent[,1])


	return( list(elaResult=elasummary, elaPlunning=redela, stability=stability ))
                    
    ## ================================ ##
    end <- proc.time()[3]
    cat(sprintf('Stability calculation done.\nElapsed time %.2f sec\n\n', end-start))                   
}
