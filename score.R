score1 <- function(coll.model, Y, E, GSP, verbose, collapse){
	model=coll.model$G
	if (verbose)
		print(model)
	if (is.vector(Y)){
		Y = matrix(Y, ncol = 1)
	}
	scores <- RMSD(model, Y, E, GSP, verbose, collapse) 
	list(RMSD=scores$RMSD, beta=scores$beta, model=scores$model)
}

score <- function(coll.models, Y, E, GSP, verbose=TRUE, collapse=TRUE, graphClass="graphNEL") {
  s=NULL
  beta=as.list(c())
  RMSDs=NULL
  for (model in coll.models){
    result = score1(model, Y, E, GSP, verbose, collapse)
    s=c(s,result$RMSD)
    beta=c(beta,list(result$beta))
    RMSDs=c(RMSDs,result$RMSD)
  }
  RMSD.diff = NULL
	if(length(s) > 1){		
		RMSD.sorted = sort(s, decreasing=TRUE)
		RMSD.diff = RMSD.sorted[[1]] - RMSD.sorted[[2]]
		if(verbose){
			cat("(Marginal RMSD difference of best vs. second best model:", RMSD.diff,")\n")
		}
  }  
  	
 	# winning model
  best <- which.min(s)
  winner <- coll.models[[best]]$G
  diag(winner) <- 1 
  beta.winner <- beta[[best]]
  RMSDs.winner <- s[best]
  #winning models
  allbest <- which(s<min(s)+1e-13)
  mods=vector("list", length(allbest)) 
  betas=vector("list", length(allbest)) 
  z=0
  if(length(s)>0){
    for(i in 1:length(allbest)){
      m=allbest[i]
      betas[[i]]=as.list(beta[[best]])
      winnerall <- coll.models[[m]]$G
      diag(winnerall) <- 1 
      mods[[i]]=winnerall
      if(sum(winnerall)>z){
        winner=winnerall
        z=sum(winnerall)
      }
      
    }
  }
  #rownames(beta.winner)=colnames(winner)
  if(graphClass == "graphNEL"){
	 	gR <- new("graphAM",adjMat=winner,edgemode="directed")  
  	gR <- as(gR,"graphNEL")    
  }else
  	gR = winner
	res <- list(graph = gR, RMSD=RMSDs.winner, beta = beta.winner, RMSD.diff = RMSD.diff, RMSDs = RMSDs, mat=winner, graphcol=collapse.cycles(winner)$G, allmat=mods, allbeta=betas)
	return(res)  
}
 