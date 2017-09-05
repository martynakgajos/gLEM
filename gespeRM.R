gespeRM <- function(Y, E, GSP, Gs=NULL, inference="greedy", verbose=FALSE, collapse = TRUE){

#------------------------------
# Sanity checks               
  if (!(inference %in% c("exhaustive", "greedy","search") )) 
    stop("\ngespeRM> argument 'inference' is not valid\n")

#------------------------------
# GREEDY                    
  if(inference == "greedy"){
    result <- gespeRM.greedy(Y, E, GSP, initial=Gs[[1]], verbose=verbose, collapse=collapse)
  }
#------------------------------
# EXHAUSTIVE SEARCH                       
  else if (inference == "exhaustive"){ 
    GEs <- enumerateGEs(E, verbose=verbose, collapse=FALSE)		
    result <- score(GEs, Y, E, GSP, verbose, collapse)	
  }

#------------------------------
# Score all given models
  else if (inference == "search"){
    if (is.null(Gs)){
    		GEs <- enumerateGEs(E, verbose=verbose, collapse=collapse)
    }else{
    		GEs = makeGEs(Gs, E, verbose= verbose, collapse=collapse)
    }
    if (length(GEs)>0){
    	result <- score(GEs, Y, E, GSP, verbose, collapse)
    }else{
    	result <- list(graph = NULL, RMSD=+Inf, beta = NULL, RMSD.diff = NULL, RMSDs = NULL)
    }	
  }

#------------------------------
# OUTPUT                       
  return(result)
}
