library(cepp)
data(Colon)

alondat <- Colon$X
alony <- Colon$Y
alonnames <- Colon$gene.names

alondatlognorm <- apply(alondat,2,log)
alondatlognorm <- apply(alondatlognorm,2,function(x) (x - mean(x))/sd(x))

zscores <- rep(NA,2000)
pvals <- rep(NA,2000)

for (i in 1:2000){
  
  X.i <- alondatlognorm[,i]
  
  my.summary <- summary(glm(y ~ X.i + 0,
                            data = data.frame(X.i,y = alony),
                            family = "binomial"))$coefficients
  
  zscores[i] <- my.summary[3]
  pvals[i] <- my.summary[4]
  
}

alonX <- alondatlognorm[,order(pvals)[1:500]]
alonnames.pared <- alonnames[order(pvals)[1:500]]

aloncass <- stan("cass_logistic.stan",
                 data = list(N = 62,
                             ncov = 500,
                             y = alony,
                             x = alonX,
                             sigma_indic = 10,
                             mu_indic = 0,
                             tau = 5),
                 chains = 4,
                 iter = 1000,
                 control = list(adapt_delta = 0.99))

## LOOCV

indices.pos <- which(alony == 1)
indices.neg <- which(alony == 0)

to.remove <- rep(NA,62)

for (i in 1:62){
  if (alony[i] == 1){
    to.remove[i] <- sample(indices.neg,1)
  } else{
    to.remove[i] <- sample(indices.pos,1)
  }
}


loolist <- list()

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

for (i in 1:62){
  
  newX <- alonX[-c(i,to.remove[i]),]
  newy <- alony[-c(i,to.remove[i])]
  
  loolist[[i]] <- stan("cass_logistic.stan",
                       data = list(N = 60,
                                   ncov = 500,
                                   y = newy,
                                   x = newX,
                                   sigma_indic = 10,
                                   mu_indic = -2,
                                   tau = 5),
                       chains = 1,
                       iter = 500,
                       control = list(adapt_delta = 0.99))
  
}

## posterior predictions

inv_logit <- function(x) 1/(1 + exp(-x))
mean.pred <- rep(NA,62)
pred <- rep(NA,250)

for (i in 1:62){
  
  beta.matrix <- extract(loolist[[i]],pars = "beta")$beta
  
  X.loo <- alonX[i,]
  
  for (j in 1:250){
    pred[j] <- rbinom(1,1,inv_logit(X.loo%*%beta.matrix[j,]))
  }
  
  mean.pred[i] <- mean(pred)
}



## loo cv for random forest

library(randomForest)

mean.pred.rf <- rep(NA,62)

for (i in 1:62){
  
  newX <- alonX[-c(i,to.remove[i]),]
  newy <- alony[-c(i,to.remove[i])]
  
  mean.pred.rf[i] <- predict(randomForest(newX,as.factor(newy),ntree = 1000),newdata = alonX[i,],type = "prob")[2]
  
}

## loo cv for lasso

library(glmnet)

mean.pred.lasso <- rep(NA,62)

for (i in 1:62){
  
  newX <- alonX[-c(i,to.remove[i]),]
  newy <- alony[-c(i,to.remove[i])]
  
  mean.pred.lasso[i] <- predict(cv.glmnet(newX,newy,family = "binomial"),newx = t(as.matrix(alonX[i,])),type = "response")
  
}

## loo cv for NN
library(MicrosoftML)

mean.pred.nn <- rep(NA,62)

nnform <- paste0("newy ~ ", paste0("X",1:500,collapse = "+"))

for (i in 1:62){
  newX <- alonX[-c(i,to.remove[i]),]
  newy <- alony[-c(i,to.remove[i])]
  
  nnfit <- rxNeuralNet(formula = nnform,data = data.frame(newX,newy),numHiddenNodes = 200,numIterations = 10000,reportProgress = 0, verbose = 0)
  
  mean.pred.nn[i] <- as.numeric(rxPredict(nnfit,data = data.frame(t(alonX[i,]),newy = alony[i])))[3]
}


## compute AUCS
library(pROC)

alon.auc.cass <- auc(roc(alony,mean.pred))
alon.auc.lasso <- auc(roc(alony,mean.pred.lasso))
alon.auc.rf <- auc(roc(alony,mean.pred.rf))
alon.auc.nn <- auc(roc(alony,mean.pred.nn))
