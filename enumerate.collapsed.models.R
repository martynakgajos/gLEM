collapse.cycles<- function( G, E=NULL, beta=NULL, GSP=NULL){
	cycle_found = TRUE
	diag(G) = 1
	if (!is.null(beta)){
	  beta=as.matrix(beta)
	  rownames(beta)=rownames(G)
	}
	if (!is.null(GSP)){
	  GSP=as.matrix(GSP)
	  rownames(GSP)=rownames(G)
	}
	while( cycle_found ){
		cycle_found=FALSE
		if ( ncol(G) > 1 ){
			for ( i in  (1: (ncol(G)-1) )){
				for ( j in ( (i+1) : ncol(G) )){
					if ( (G[i,j]==1) & (G[j,i]==1) ){
						#collapse
						node_i_name = colnames(G)[i]
						node_j_name = colnames(G)[j]
						new_node_name = paste( c( node_i_name, "_", node_j_name), collapse="")
						
						in.edges = as.numeric( ( G[, node_i_name] + G[, node_j_name]) > 0 )
						out.edges = G[i, ]
						# add new collapsed node
						G = cbind(in.edges, G)
						G = rbind(c(1, out.edges), G)
						colnames(G) = rownames(G) = c(new_node_name, colnames(G)[2:ncol(G)]) 

						## add the node to E 
						if (! is.null(E)){
							in.experiments <- as.numeric( ( E[, node_i_name] + E[, node_j_name]) > 0 )
							E = cbind(in.experiments, E)
							colnames(E) <- c(new_node_name, colnames(E)[2:ncol(E)])
						}
						
						## make the total contribution
						if (!is.null(beta)){
						  #beta = rbind(sum(beta[i,], beta[j,]), beta )
							beta = rbind(sum(beta[i], beta[j]), beta )
							rownames(beta ) = c(new_node_name, rownames( beta)[2:nrow(beta)] )
						}
						if (!is.null(GSP)){
						  GSP = rbind(sum(GSP[i], GSP[j])/2., GSP )
						  rownames(GSP) = c(new_node_name, rownames(GSP)[2:nrow(GSP)] )
						  GSP = GSP[-( which(rownames(GSP) == node_i_name)), ,drop=F]
						  GSP = GSP[-( which(rownames(GSP) == node_j_name)), ,drop=F]
						}

						# remove the nodes i and j
						G = G[ -( which(rownames(G) == node_i_name)), , drop=F ]
						G = G[ -( which(rownames(G) == node_j_name)), , drop=F ]
						G = G[ , -( which(colnames(G) == node_j_name)), drop=F ]
						G = G[ , -( which(colnames(G) == node_i_name)), drop=F ]

						if (! is.null(E)){
							E = E[,-( which(colnames(E) == node_i_name)), drop=F ]
							E = E[,-( which(colnames(E) == node_j_name)), drop=F ]
						}	
						
						if (!is.null(beta)){
							beta = beta[-( which(rownames(beta) == node_i_name)), ,drop=F]
							beta = beta[-( which(rownames(beta) == node_j_name)), ,drop=F]
						}
	
						cycle_found = TRUE
						break();
					} 					
				}
				if (cycle_found) break()
			}
		}
	}
	res = list(G = G, E = E, beta = beta, GSP=GSP)

	res
}

makeGEs <- function(Gs, E, verbose=FALSE, collapse=TRUE){
	### Remove the model which will be collapsed into 1 node
	alln1 <- sapply( Gs, function(G){ !all(G==1) })
	Gs <- Gs[alln1]
	if (collapse){
		#Gs	<- lapply(Gs, collapse.cycles, E = E )
	  Gs	<- lapply(Gs, collapse.cycles)
	}else{
		Gs = lapply( Gs, function(G){ list(G=G, E=E) })
	}
	Gs
}

enumerateGEs<-function(E, verbose=FALSE, collapse=TRUE){
	Gs<-enumerate.models(ncol(E), colnames(E), trans.close=TRUE, verbose=verbose)
	makeGEs(Gs=Gs, E=E, verbose=verbose, collapse=collapse)
}

