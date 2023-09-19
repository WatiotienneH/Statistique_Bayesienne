#Génération d'un jeu de données
x = rnorm(50)
y = 2 + 3*x + rnorm(50)


#Approche fréquentiste
REG=lm(y~x)
summary(REG)
plot(x,y)
curve(REG$coefficients[1]+REG$coefficients[2]*x,add = T)
confint(REG)


#approche bayésienne
library(Bolstad)
bayes.lin.reg(y,x) #par défaut les priors sont non informatifs


##############################################################
#loi a priori informative

library(wooldridge)
data = wooldridge::wine
data;
modele=lm(data$heart~data$alcohol)
summary(modele)
plot(data$alcohol,data$heart)
curve(239.147-19.683*x,add=T)




bayesien_modele=bayes.lin.reg(
  data$heart,
  data$alcohol,
  slope.prior = "normal",
  intcpt.prior = "flat",
  mb0=-15,
  sb0=5,
  ma0=0,
  sa0=0,
  sigma=NULL,
  alpha=0.05,
  plot.data=TRUE,
  pred.x = NULL,
  
)
