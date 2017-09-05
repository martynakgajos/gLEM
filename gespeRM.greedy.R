get.insertions = function(Phi){
    idx = which(Phi == 0)
    models = list()
    if(length(idx) > 0){
        for(i in 1:length(idx)){ # test all possible new edges
            Phinew = Phi
            Phinew[idx[i]] = 1
            Phinew = transitive.closure(Phinew, mat=TRUE,loops=TRUE) 
            models[[i]] <- Phinew
        }
    } 
    models       
}

delete.edge.transclose <- function( Phinew,  rowno , colno){
	Phinew[rowno, colno] = 0
	### update that the parents of gene rowno also remove edges to gene colno
	# get the parents of rowno
	parents = which( Phi[,rowno] == 1 )
	# remove edges parents -- colno
	Phinew[parents, colno] = 0
	# still some parents could feed into colno via independent routes. For them we need to restore their edges to colno by transitive closure
	Phinew= transitive.closure(Phinew, mat=TRUE,loops=TRUE) 
	Phinew	
}

get.deletions = function(Phi){
    Phi = Phi - diag(ncol(Phi))
    idx = which(Phi == 1)
    models = list()
    if(length(idx) > 0){
        for(i in 1:length(idx)){ # test all possible edge deletions
            Phinew = Phi
            ### retreive the column at which we made the switch to 0
            colno = ceil(idx[i]/nrow(Phi))
            ### retreive the row at which we made the switch to 0. This means that the parent rowno no longer has an edge to the gene colno
			      rowno = idx[i] - (colno-1)*nrow(Phi)
			
			if (colno!=rowno){ # we do not make an attempt to remove the diagonal
				Phinew = delete.edge.transclose( Phinew,  rowno , colno)
				diag(Phinew) = 1
				models[[i]] <- Phinew
			}
      }
    } 
    models       
}

get.reversions = function(Phi){
    idx = which(Phi + t(Phi) == 1, arr.ind=TRUE)
    models = list()
    if(NROW(idx) > 0){
        for(i in 1:NROW(idx)){ # test all possible edge reversions
            Phinew = Phi
            Phinew[idx[i,1],idx[i,2]] = 0
            Phinew[idx[i,2],idx[i,1]] = 1
            diag(Phinew) = 1
            models[[i]] <- Phinew
        }
    } 
    models       
}

gespeRM.greedy <- function(Y, E, GSP, initial=NULL, verbose=TRUE, collapse=TRUE){ 
    n <- length(GSP) 
    if (verbose){
    	cat("Greedy hillclimber for",n,"genes...\n\n")
    }
    if (is.null(initial)){
    	Phi <- matrix(0,nrow=n,ncol=n)
    	colnames(Phi)=colnames(E)
    }else{
    	Phi = initial
    }
    diag(Phi) <- 1    
    sco0 <- gespeRM(Y, E, GSP, Gs=list(Phi), inference = "search", verbose=verbose, collapse=FALSE)$RMSD
    finished <- FALSE
    i = 1
    while(!finished){
        models <- list()
#       propose new edges     
        models = get.insertions(Phi)
        models <- unique(models)
        if(verbose)
            cat(length(models), " local models to test ...\n")
        if(length(models) > 0){
            sconew <- gespeRM(Y, E, GSP, Gs=models, inference="search", verbose=verbose, collapse=FALSE)          
            if(min(sconew$RMSD) <= sco0){
                if(verbose)
                    cat("step",i,"--> Edge added(!), removed or reversed\n")
                sco0 <- min(sconew$RMSD)
                Phi <- as(sconew$graph,"matrix")            
            }
            else # otherwise no improving edge could be inserted
                finished <- TRUE
        }else
            finished <- TRUE    
    	i = i+1
    }
    if (verbose){
    	cat("lem greedy made",i,"iterations\n\n")
    }
    res <- gespeRM(Y, E, GSP, Gs=list(Phi), inference="search", verbose=verbose, collapse=FALSE)
    
    if(verbose)
        cat("RMSD of the model = ",res$RMSD,"\n\n")
    return(res)
}
