#chargement des bibliothèques
library(splines)
library(fda)
library(MASS)

#On génére un jeu de données
n <- 100
x <- runif(n, 0, 1)
y <- sin(2 * pi * x) + rnorm(n, 0, 0.2)# On ajoute un bruit
plot(x,y) #on visualise nos données
curve(sin(2 * pi * x),add=T,col="green") #la vraie fonction

#Création de la matrice B et Omega
r=4 #spline cubique
M=(n-2)+r #nombre de fonctions de bases
z=x[order(x)]
bspl=create.bspline.basis(range(z),nbasis=M,norder=r,z) #création d’une base BSpline - range(z) = intervalle sur lequel les fonctions doivent être évaluées 
Omega=bsplinepen(bspl,Lfdobj=2,rng=bspl$rangeval)  #création de la matrice de pénalisation (pénalité de lissage pour fonctions exprimées en base B-spline)
B=bsplineS(x,z,nord=r)


#paramétre de lissage 
lambda=0.0001 #déterminer empiriquement
S=B%*%solve((t(B)%*%B+lambda*Omega))%*%t(B) #la matrice S

sigma2=0.2 #On suppose sigma2 connue


f=mvrnorm(2000,S%*%y,sigma2*S) #On génére 2000 échantillons de f|y
f_med=apply(f,2,median) #on récupére la médiane
f_lower=apply(f,2,quantile,probs=0.025)
f_upper=apply(f,2,quantile,probs=0.975)
plot(x,y)
lines(z,f_med[order(x)])
lines(z,f_lower[order(x)],lty=3,col='green')
lines(z,f_upper[order(x)],lty=3,col='green')
curve(sin(2 * pi * x),add=T,col="red")
legend("topright", legend = c("f*(x)","intervalle de prédiction 95%"), col = c("black","green"), lty = c(1,3))  
legend("bottomright", legend = "f(x)", col = "red", lty = 1)  

