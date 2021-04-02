##Données##

N <- 27  
x <- c(1, 1.5, 1.5, 1.5, 2.5, 4, 5, 5, 7, 8, 8.5, 9, 9.5, 9.5, 10, 
       12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29, 31.5)
Y <- c(1.8, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 
       2.26, 2.4, 2.39, 2.41, 2.5, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 
       2.47, 2.64, 2.56, 2.7, 2.72, 2.57)

##Initialisation et paramètres des lois a priori##

alpha <- 2.5
beta <- 1
tau <- 1
gamma <- 0.9


s <- 10^4
a <- b <- 0.001
sd.prop <- c(0.07,0.5,1,0.1)

Nchain <- 100000
chain <- matrix(NA, Nchain+1, 4)


##Création de la chaine vide et du taux d'acceptation##

chain[1,] <- c(alpha, beta, tau, gamma)
acc.rates <- c(0,0,0,0)

for (iter in 1:Nchain){
  
  current <- chain[iter,]
  
  ##Nouveau alpha##
  
  prop.alpha <- rlnorm(1, meanlog = log(current[1]), sd.prop[1])
  ratio <- prop.alpha/current[1]
  acc.top <- sum(dnorm(Y,prop.alpha - current[2]*current[4]^x, 1/sqrt(current[3]), log=TRUE)) +
                   dnorm(prop.alpha, 0, s, log=TRUE)
  acc.bottom <- sum(dnorm(Y,current[1] - current[2]*current[4]^x, 1/sqrt(current[3]), log=TRUE)) +
                      dnorm(current[1], 0, s, log=TRUE)
  
  acc <- exp(acc.top - acc.bottom)*ratio
  
  if (runif(1)<acc){
    current[1] <- prop.alpha
    acc.rates[1] <- acc.rates[1] + 1
  }
  
  
  
  ##Nouveau beta##
  
  prop.beta <- rlnorm(1, meanlog = log(current[2]), sd.prop[2])
  ratio <- prop.beta/current[2]
  acc.top <- sum(dnorm(Y,current[1] - prop.beta*current[4]^x, 1/sqrt(current[3]), log=TRUE)) +
                   dnorm(prop.beta, 0, s, log=TRUE)
  acc.bottom <- sum(dnorm(Y,current[1] - current[2]*current[4]^x, 1/sqrt(current[3]), log=TRUE)) +
                      dnorm(current[2], 0, s, log=TRUE)

  acc <- exp(acc.top - acc.bottom)*ratio

  if (runif(1)<acc){
    current[2] <- prop.beta
    acc.rates[2] <- acc.rates[2] + 1
  }
  
  
  ##Nouveau tau##
  
  prop.tau <- rlnorm(1, meanlog = log(current[3]), sd.prop[3])
  ratio <- prop.tau/current[3]
  acc.top <- sum(dnorm(Y,current[1] - current[2]*current[4]^x, 1/sqrt(prop.tau), log=TRUE)) +
                   dgamma(prop.tau, a, b, log=TRUE)
  acc.bottom <- sum(dnorm(Y,current[1] - current[2]*current[4]^x, 1/sqrt(current[3]), log=TRUE)) +
                      dgamma(current[3], a, b, log=TRUE)
  
  acc <- exp(acc.top - acc.bottom)*ratio
  
  if (runif(1)<acc){
    current[3] <- prop.tau
    acc.rates[3] <- acc.rates[3] + 1
  }
  
  
  ##Nouveau gamma##
  
  
  prop.gamma <- rlnorm(1, meanlog = log(current[4]), sd.prop[4])
  ratio <- prop.gamma/current[4]
  acc.top <- sum(dnorm(Y,current[1]-current[2]*prop.gamma^x, 1/sqrt(current[3]), log=TRUE)) + 
    dunif(prop.gamma, 0.5, 1, log=TRUE)
  acc.bottom <- sum(dnorm(Y,current[1]-current[2]*current[4]^x, 1/sqrt(current[3]), log=TRUE)) + 
    dunif(current[4], 0.5, 1, log=TRUE)
  
  acc <- exp(acc.top - acc.bottom)*ratio
  
  if ((runif(1)<acc)&(prop.gamma<1)){
    current[4] <- prop.gamma
    acc.rates[4] <- acc.rates[4] + 1
  }
  chain[iter+1,] <- current
}


## Affichage des observations##
plot(x,Y,main='Poids et âge des 27 dugongs', xlab='Age (en années)', ylab='Taille (en m)')


## taux d'acceptation##
acc.rates <- acc.rates/Nchain
acc.rates 


##Chaine de alpha avec ACF et tronquage##
plot(1500:length(chain[,1]), chain[-(1:1500),1], main='Chaine de alpha', xlab='Itérations', ylab="alpha")
acf(chain[-(1:1500),1], main='ACF de alpha')

##Tronquage des chaines##
chaintronq.alpha <- chain[-(1:1500),1][seq(1,length(chain[-(1:1500),1]), by=10)]
chaintronq.beta <- chain[-(1:1500),2][seq(1,length(chain[-(1:1500),1]), by=10)]
chaintronq.tau <- chain[-(1:1500),3][seq(1,length(chain[-(1:1500),1]), by=10)]
chaintronq.sigma <- 1/sqrt(chaintronq.tau)
chaintronq.gamma <- chain[-(1:1500),4][seq(1,length(chain[-(1:1500),1]), by=10)]
plot(1:length(chaintronq.alpha), chaintronq.alpha, main='Chaine de alpha (un élément sur 10)', xlab='Itérations', ylab="alpha")
acf(chaintronq.alpha, main='ACF de alpha (chaine tronquée)')

##Histogrammes##
hist(chaintronq.alpha, breaks=50, freq=FALSE, main="Densité de alpha", xlab='alpha', col='orange')
densite <- density(chaintronq.alpha) 
lines(densite, col = "red",lwd=3)

hist(chaintronq.beta, breaks=50, freq=FALSE, main="Densité de beta", xlab='beta', col='orange')
densite <- density(chaintronq.beta) 
lines(densite, col = "red",lwd=3)

hist(chaintronq.sigma, breaks=50, freq=FALSE, main="Densité de sigma", xlab='sigma', col='orange')
densite <- density(chaintronq.sigma) 
lines(densite, col = "red",lwd=3)

hist(chaintronq.gamma, breaks=50, freq=FALSE, main="Densité de gamma", xlab='gamma', col='orange')
densite <- density(chaintronq.gamma) 
lines(densite, col = "red",lwd=3)


##Données statistiques##
mean(chaintronq.alpha)
sd(chaintronq.alpha)
median(chaintronq.alpha)
quantile(chaintronq.alpha, 0.025)
quantile(chaintronq.alpha, 0.975)

mean(chaintronq.beta)
sd(chaintronq.beta)
median(chaintronq.beta)
quantile(chaintronq.beta, 0.025)
quantile(chaintronq.beta, 0.975)

mean(chaintronq.sigma)
sd(chaintronq.sigma)
median(chaintronq.sigma)
quantile(chaintronq.sigma, 0.025)
quantile(chaintronq.sigma, 0.975)

mean(chaintronq.gamma)
sd(chaintronq.gamma)
median(chaintronq.gamma)
quantile(chaintronq.gamma, 0.025)
quantile(chaintronq.gamma, 0.975)

alpha <- mean(chain[-(1:1500),1])
beta <- mean(chain[-(1:1500),2])
gamma <-  mean(chain[-(1:1500),4])


## Observations et courbe approximée##

plot(x,Y,main='Poids et âge des 27 dugongs', xlab='Age (en années)', ylab='Taille (en m)')
lines(seq(0,35,0.1), alpha - beta*gamma^seq(0,35,0.1), col='red', lwd=3)

##Calcul de sigma_tild##

alpha.mean <- mean(chaintronq.alpha)
beta.mean <- mean(chaintronq.beta)
sigma.mean <- mean(chaintronq.sigma)
gamma.mean <- mean(chaintronq.gamma)

sigma.tild = sqrt(mean((Y-alpha.mean+beta.mean*gamma.mean^x)^2))
sigma.tild
sigma.mean