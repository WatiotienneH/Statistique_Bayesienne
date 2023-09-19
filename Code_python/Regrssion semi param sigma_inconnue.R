#chargement des bibliothèques
library(splines)
library(fda)
library(MASS)
library(invgamma)

#On génére un jeu de données pour verifier notre code avant de l'appliquer sur un vrai jeu de données

n <- 100
x <- runif(n, 0, 1)
y <- sin(2 * pi * x) + rnorm(n, 0, 0.2) # On ajoute un bruit
plot(x,y) #on visualise nos données
curve(sin(2 * pi * x),add=T,col="green") #la vraie fonction


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


tau2=8 #déterminer empiriquement




N=1000 # Nombre d'itération

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




#On récupére la derniére valeur de sigma2 pour construire un f
S=B%*%solve((t(B)%*%B+lambda[N]*Omega))%*%t(B)
f_final=apply(f_iter[floor(0.8*N):N,],2,median)
f_lower=apply(f_iter[floor(0.8*N):N,],2,quantile,probs=0.025)
f_upper=apply(f_iter[floor(0.8*N):N,],2,quantile,probs=0.975)
plot(x,y)
lines(z,f_final[order(x)])
lines(z,f_lower[order(x)],lty=1,col='blue')
lines(z,f_upper[order(x)],lty=1,col='blue')
curve(sin(2 * pi * x),add=T,col="red")
legend("topright", legend = c("f*(x)","intervalle de prédiction 95%"), col = c("black","blue"), lty = c(1,1))  
legend("bottomright", legend = "f(x)", col = "red", lty = 1)  


#représentation des chaines
plot(1:N,sigma2)
