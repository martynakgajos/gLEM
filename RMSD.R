RMSD <- function(G, Y, E, GSP, verbose=FALSE, collapse=TRUE) {
  require(matlib)
  
  #stworzenie macierzy S
  S=E
  for (i in 1:dim(S)[1]){
    buf=which(S[i,]==1)
    for (j in buf){
      S[i,which(G[j,]==1)]=1
    }
  }
  
  GSP=as.matrix(GSP)
  rownames(GSP)=rownames(G)
  if (collapse==TRUE|det(G)==0){
    res=collapse.cycles(G, S, GSP=GSP)
    GSP=res$GSP
    S=res$E
    G=res$G
  }
  #lm.model=lm(GSP~0+G)
  #beta=lm.model$coefficients
  #beta[which(is.na(beta))] <- 1
  beta=solve(G,GSP)
  Yprim=S %*% beta
  RMSD=sum(sqrt((Yprim-Y)**2))
  res=list(RMSD=RMSD, beta=beta,model=G)
  # if (RMSD<9.8){
  #   gR <- new("graphAM",adjMat=G,edgemode="directed")  
  #   gR <- as(gR,"graphNEL")
  #   B$graph=gR
  #   png(filename=paste(runif(1,0,100),'.png'))
  #   plot.nem(B)
  #   dev.off()
  # }
  # print("G")
  # print(G)
  # print("S")
  # print(S)
  # print("GSP")
  # print(GSP)
  # print("beta")
  # print(beta)
  # print("RMSD")
  # print(RMSD)

  res
}
