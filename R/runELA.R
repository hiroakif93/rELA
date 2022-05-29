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

	speciesName <- rownames(J)
	
	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Stable state estimatin	
	minsets = SSestimate(alpha, J, itr = SS.itr)
	minsets = unique(round(minsets, digits = 8))
	minsets <- minsets[order(minsets[, ncol(minsets)]), ]
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
	
	FindTip <- c()
	for (k in 1:nrow(comb)) {
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
    #for.parallel(thread)            	
    start <- proc.time()[3]

	speciesName <- rownames(J)
	
	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Stable state estimatin	
	minsets = SSestimate(alpha, J, itr = SS.itr)
	minsets = unique(round(minsets, digits = 8))
	minsets <- minsets[order(minsets[, ncol(minsets)]), ]
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
ELplunning <- function(elasummary, th=0.9){
    
	stateInfo <- elasummary[[1]]
	elaMat <- elasummary[[2]]
	
	## ========================================== ##
	## -- Make incidence matrix
	
	IM <- matrix(Inf, nrow=sum(stateInfo$state=='Stable state'), ncol=nrow(elaMat)+1,
	             dimnames=list(stateInfo$stateID[which(stateInfo$state=='Stable state')]))
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
	
	return(list(matrix=elaMattmp, log=log[,c(1,ncol(log))]))
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
						  seitr=1000, convTime=10000){
						  	
    start <- proc.time()[3]
    ## ================================ ##
    sampleSS <- t(apply(data, 1, SteepestDescent_cpp, alpha=alpha, beta=J))
    ssid <- apply(sampleSS[,-ncol(sampleSS)], 1, 
                  function(x){ paste0(x, collapse='') })
    ssenergy <- sampleSS[, ncol(sampleSS)]
    energy <- apply(ocdata, 1, cEnergy, alpha= alpha, beta=J)
    
    energy.gap <- energy-ssenergy
    ssent <- SSentropy_cpp(uoc=ocdata, ss=unique(sampleSS[,-ncol(sampleSS)]),
          				alpha= alpha, beta=J, 
          				seitr=seitr, convTime=convTime)
          	  
    return( data.frame(StableState.id=ssid, SSenergy=ssenergy,
                       energy.gap=energy.gap, SSentropy= ssent[,1]))
                       
    ## ================================ ##
    end <- proc.time()[3]
    cat(sprintf('Stability calculation done.\nElapsed time %.2f sec\n\n', end-start))                   
}
