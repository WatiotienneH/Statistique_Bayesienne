library(splines)
library(fda)
library(MASS)
library(invgamma)

chemin_fichier <-"C:/Users/hrywa/Desktop/Donnees_temperature/TX_STAID000037.txt"
donnees <- read.table(chemin_fichier, sep = ",",header=TRUE ,skip = 20)

#on retire les valeurs manquantes
donnees_filtrees <- subset(donnees, TX != -9999)
# Indices des valeurs conservées sur trois
indices_valeurs_conservees <- seq(1, length(donnees_filtrees$TX), by = 40)
y<- 0.1*donnees_filtrees$TX[seq(1, length(donnees_filtrees$TX), by = 40)]
n=length(y)
x=1:n
plot(x,y)

#Construction de B et Omega
r=4  #spline cubique
M=(n-2)+r #nombre de fonctions de bases
z=x[order(x)]
bspl=create.bspline.basis(range(z),nbasis=M,norder=r,z) #création d’une base BSpline
Omega=bsplinepen(bspl,Lfdobj=2,rng=bspl$rangeval)  #création de la matrice de pénalisation (pénalité de lissage pour fonctions exprimées en base B-spline)
B=bsplineS(x,z,nord=r)



#hyperparamétre de la loi a priori sur sigma^2 suit une invgamma(a,b)
a=0.01
b=0.01
curve(dinvgamma(x,a,b),xlim = c(0,10),main="Loi a priori de sigma^2")


tau2=1 #déterminer empiriquement




N=50 # Nombres d'itérations

#initialisation

f_iter=matrix(y,N,n) #on stocke les vecteurs f dans une matrice 
sigma2=rep(0.2,N)
lambda=rep(0.0015,N)
# Algorithme Gibbs

for(i in 2:N){
  S=B%*%solve((t(B)%*%B+lambda[i-1]*Omega))%*%t(B)#on calcul S(lambda)
  f=mvrnorm(1,S%*%y,sigma2[i-1]*S)  #On simule f|sigma^2,Y en utilisant le précédent sigma^2
  f_iter[i,]=f                  #On stocke le f dans la matrice
  sigma2_value=rinvgamma(1,n/2+a,(1/2)*norm(y-f_iter[i,],type="2")^2+b) #on simule sigma2|Y,f en utilisant le f que l'on vient de simuler
  sigma2[i]=sigma2_value      #on stocke la valeur mediane
  lambda[i]=sigma2[i]/tau2
}


plot(x, y, pch = ".", col = "transparent")
f_final=apply(f_iter[floor(0.8*N):N,],2,median)
f_lower=apply(f_iter[floor(0.8*N):N,],2,quantile,probs=0.025)
f_upper=apply(f_iter[floor(0.8*N):N,],2,quantile,probs=0.975)
lines(z,f_final[order(x)])
lines(z,f_lower[order(x)],lty=1,col='blue')
lines(z,f_upper[order(x)],lty=1,col='blue')
legend("topright", legend = c("f*(x)","intervalle de prédiction 95%"), col = c("black","blue"), lty = c(1,1))  


