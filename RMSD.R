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
  #obliczenie bet
  beta=solve(G,GSP)
  Yprim=S %*% beta
  #wyznaczenie RMSD
  RMSD=sum(sqrt((Yprim-Y)**2))
  res=list(RMSD=RMSD, beta=beta,model=G)
  
  res
}
