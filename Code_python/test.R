library(invgamma)
library(ggplot2)
# Modele bayesien

# Xi | mu, sigma^2  ~ N( mu, sigma^2)
#    mu|sigma^2 ~ N(a, sigma^2 /b)
#       sigma^2 ~ IG(c,d)





# les hyperparamétres
a <- 8       
b <- 0.5
c <- 1
d <- 1


par(mfrow=c(3,2))
chaine <- matrix(c(1,1),nc=2,nr=1) #initialisation 
for(j in 1:3){
  N <- c(10,100,1000)
  n=N[j]
  X=rnorm(n,9,1.5)
  Xn=mean(X)
  sum_X=sum(X)
  SEM=sum((X -mean(X))^2)
  chaine <- matrix(c(1,1),nc=2,nr=1)
  
  for (i in 2:n){
    simga_2 <- rinvgamma(1,c+n/2,d+0.5*(SEM+(n*b)*(Xn-a)^2/(2*(n+b))))
    mu <- rnorm(1,(sum_X+a*b)/(n+b),simga_2/(a+b))
    chaine <- rbind(chaine,c(simga_2,mu))
  }
  hist(chaine[2:n,1],main=paste("Distribution de la loi a posteriori de sigma^2 n=",N[j]),freq=F,col="grey",xlim=c(0,3),xlab=" ")
  curve(dinvgamma(x,c,d),col="red",add=T)
  abline(v=SEM/n,col="blue")
  legend("topleft",legend = c("prior","EMV"), col = c("red","blue"), lty = c(1,1))
  hist(chaine[2:n,2],freq = F,main=paste("Distribution de la loi a posteriori de mu n=",N[j]),col="grey",xlim=c(5,15),xlab=" ")
  curve(dnorm(x,a,1.5/b),col="red",add=T)
  abline(v=Xn,col="blue")
  legend("topleft",legend = c("prior","EMV"), col = c("red","blue"), lty = c(1,1))
}









