sampleRndGraph = function(E, collapse=TRUE){	
	
	allone=TRUE
	while (allone){
		G = sampleRndNetwork(colnames(E))
		diag(G)=1
		allone=all(G==1)
	}
	if (collapse){
		GE	<-  collapse.cycles(G=G, E = E )
	}else{
		GE <- list(G  = G, E = E)
	}
	
	GE

}
