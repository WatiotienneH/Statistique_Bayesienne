library(splines)
library(fda)
library(MASS)
library(invgamma)

chemin_Lyon <-"C:/Users/hrywa/Desktop/Donnees_temperature/TX_STAID000037.txt"
chemin_BRINDISI<-"C:/Users/hrywa/Desktop/Donnees_temperature/TX_STAID000174.txt"
chemin_AKUREYRI<-"C:/Users/hrywa/Desktop/Donnees_temperature/TX_STAID002943.txt"
chemin_HELSINKI<-"C:/Users/hrywa/Desktop/Donnees_temperature/TX_STAID007015.txt"
chemin_STERLEGOVA<-"C:/Users/hrywa/Desktop/Donnees_temperature/TX_STAID008044.txt"


a=0.01
b=0.01
curve(dinvgamma(x,a,b),xlim = c(0,10),main="Loi a priori de sigma^2")

donnees <- read.table(chemin_Lyon, sep = ",",header=TRUE ,skip = 20)
#on retire les valeurs manquantes
donnees_filtrees <- subset(donnees, TX != -9999)
# Indices des valeurs conservées sur trois
# Nombre de groupes de 30 valeurs
nb_groupes <- length(donnees_filtrees$TX) %/% 365


# Vecteur pour stocker les moyennes
vecteur_moyennes <- numeric(nb_groupes)

# Calculer les moyennes de chaque groupe de 30 valeurs
for (i in 1:nb_groupes) {
  debut <- (i - 1) * 365 + 1
  fin <- i * 365
  vecteur_groupe <- donnees_filtrees$TX[debut:fin]
  moyenne_groupe <- mean(vecteur_groupe)
  vecteur_moyennes[i] <- moyenne_groupe
}

y=0.1*vecteur_moyennes

n=length(y)
x=1951:(n-1+1951)
plot(x,y)


#Construction de B et Omega
r=4  #spline cubique
M=(n-2)+r #nombre de fonctions de bases
z=x[order(x)]
bspl=create.bspline.basis(range(z),nbasis=M,norder=r,z) #création d’une base BSpline
Omega=bsplinepen(bspl,Lfdobj=2,rng=bspl$rangeval)  #création de la matrice de pénalisation (pénalité de lissage pour fonctions exprimées en base B-spline)
B=bsplineS(x,z,nord=r)

tau2=0.009 #déterminer empiriquement

N=5000 # Nombres d'itérations

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



f_final=apply(f_iter[(0.8*N):N,],2,median)
f_lower=apply(f_iter[(0.8*N):N,],2,quantile,probs=0.025)
f_upper=apply(f_iter[(0.8*N):N,],2,quantile,probs=0.975)
plot(x, y,ylab = "Moyenne des température max",xlab="Années")
lines(z,f_final[order(x)])
lines(z,f_lower[order(x)],lty=1,col='blue')
lines(z,f_upper[order(x)],lty=1,col='blue')
legend("topright", legend = c("f*(x)","intervalle de prédiction 95%"), col = c("black","blue"), lty = c(1,1))  

