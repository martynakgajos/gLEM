library(matlib)
library(nem)
library(lem)

simulate.exp <- function( n, rep.no ){

  m = (n + (n-1)*n/2 )*rep.no
  E = matrix(0, ncol = n, nrow= m)
  colnames(E) = letters[1:n]
  e.double = rep.no*n + 1
  
  for( i in seq(1, n-1)){
    for (j in seq(i+1, n)){
      for (r in seq(1, rep.no)){
        E[e.double, i ] = 1
        E[e.double, j ] = 1
        e.double = e.double+1
      }
    }
  }
  
  e.single = 1
  for( i in seq(1, n) ){
    for ( r.s in seq( 1, rep.no ) ){
      E[e.single, i] = 1
      e.single = e.single+1
    }
  }
  E
}

simulate.data <- function(G, E, b, sigma,Y.no=1){
  #cycless=collapse.cycles(G, E)
  #G=cycless$G
  #E=cycless$E
  #Y.no = 1
  n = ncol(E)
  m = nrow(E)
  #betas =	sapply(1:Y.no, function(i){	be = abs(rnorm(n, 0, sd = sqrt(1/b))); names(be) = colnames(E); be } )
  betas =	sapply(1:Y.no, function(i){	be = rnorm(n, 0, sd = sqrt(1/b)); names(be) = colnames(E); be } )
  epsilon = matrix( rnorm( m*Y.no, 0, sigma), ncol=Y.no)
  X = ( (E%*%G)>0)
  
  Y = X%*%betas + epsilon
  return(list(Y=Y, betas=betas, E=E, G=G))
}

solve.simulated<- function(E,Y,i){
  A=as.matrix(E[1:ncol(E),1:ncol(E)])
  B=as.vector(Y[1:ncol(E)])
  lm.model=lm(B~0+A)
  GSP=lm.model$coefficients
  res=gespeRM(Y,E,GSP,inference=i,collapse=FALSE)
  return(res)
}

b=1
sigma=0.00
n=3
rep=1
i='exhaustive'
E=simulate.exp(n,rep)
z=1
for (M in enumerate.models(n)){
#   #res=simulate.data(M,E,b,sigma)
#   #E=res$E
#   #Y=res$Y
#   #M=res$G
#   #wynik=solve.simulated(E,Y,i)
#   print(z)
#   print(M)
#   #print(wynik$allmat)
#   #print(wynik$RMSD)
#   z=z+1
# }
A=enumerate.models(n)
M=A[[19]]
res=simulate.data(M,E,b,sigma)
E=res$E
Y=res$Y
M=res$G
wynik=solve.simulated(E,Y,i)
print("M")
print(M)
print(wynik$allmat)
print(wynik$RMSD)