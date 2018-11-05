# simulation data ----

int.20 <- 0.6
int.70 <- -0.3
int.120 <- 0.4


beta.20 <- c(rep(0,5),rep(2,5),rep(0,5),(rep(-1,5) + rnorm(5,mean = 0,sd = 0.1)))

beta.70 <- c(rep(0,15),rep(2,5),rep(0,10),rep(2,5),(rep(-1,5) + rnorm(5,mean = 0,sd = 0.1)),rep(0,25),-0.5,-0.5,3,-0.5,-0.5)

beta.120 <- c(rep(0,25),rep(1,5),rep(0,25),(rep(2,5) + rnorm(5,mean = 0,sd = 0.1)), rep(0,25),-1,-1,-1,-2,-1,0.5,.5,.5,2,.5,rep(0,25))

cov.20 <- randomLHS(100,20)

cov.70 <- randomLHS(100,70)

cov.120<- randomLHS(100,120)

y.20 <- int.20 + cov.20 %*% beta.20 + rnorm(100,mean = 0,sd = 1)

y.70 <- int.70 + cov.70 %*% beta.70 +rnorm(100,mean = 0, sd = 0.5)

y.120 <- int.120 + cov.120 %*% beta.120 + rnorm(100,mean = 0, sd = 0.5)




# 


# run some fitting algorithms

#### p = 20, n = 100 ----

grplassofit <- cvSGL(data = list(x = cov.20,y = y.20),
                     index = c(rep(1,5),rep(2,5),rep(3,5),rep(4,5)),
                     maxit = 10000)

hsfit <- stan_glm(y.20 ~ . ,
                  data = data.frame(cov.20,y.20),
                  prior = hs(df = 2,global_scale = 0.2),
                  family = gaussian(),
                  prior_intercept = normal(0,10),
                  adapt_delta = 0.9999)

median.pars.hs <- apply(as.matrix(hsfit),2,median)[2:21]

lmfit <- lm(y.20 ~ . , data = data.frame(cov.20,y.20))

glmfit <- cv.glmnet(x = cov.20,y = y.20)

cassfit <- stan("simcasshered.stan",
                data = list(N = 100,
                            ncov = 20,
                            y = y.20,
                            X = cov.20,
                            sigma_indic = 10,
                            mu_indic = 0,
                            ngroup = 4),
                control = list(adapt_delta = 0.9999),
                chains = 4, cores = 4, iter = 10000)

median.pars.cass <- summary(cassfit,pars = "beta")$summary[,6]


newdf <- data.frame(actual = rep(beta.20,5),
                    fitted = c(grplassofit$fit$beta[,8],
                               median.pars.hs,
                               coef(lmfit)[-1],
                               coef(glmfit)[-1],
                               median.pars.cass))
newdf <- data.frame(melt(newdf),class = rep(rep(c("SGL","Horseshoe","OLS","LASSO","LN-CASS"),each = 20),2))

ggplot(data = newdf, 
       aes(x = rep(1:20,10),y = value,col = variable)) + 
  geom_point(size = 1) + 
  geom_line() + 
  facet_wrap(~ class,nrow = 1,scales = "free_y") + 
  scale_color_manual(labels=c("Actual", "Fitted"),
                     values = c("blue","red")) +
  labs(x = "Coefficient index", 
       y = "Value of coefficient",
       color = "")


#### p = 70, n = 100 ----

grplassofit1 <- cvSGL(data = list(x = cov.70,y = y.70),
                     index = rep(1:14,each = 5),
                     maxit = 10000)

hsfit1 <- stan_glm(y.70 ~ . ,
                  data = data.frame(cov.70,y.70),
                  prior = hs(df = 2),
                  family = gaussian(),
                  prior_intercept = normal(0,10),
                  adapt_delta = 0.9999)

median.pars.hs1 <- apply(as.matrix(hsfit1),2,median)[2:71]


lmfit1 <- lm(y.70 ~ . , data = data.frame(cov.70,y.70))

glmfit1 <- cv.glmnet(x = cov.70,y = y.70)

cassfit1 <- stan("simcasshered.stan",
                data = list(N = 100,
                            ncov = 70,
                            y = y.70,
                            X = cov.70,
                            sigma_indic = 10,
                            mu_indic = 0,
                            ngroup = 14),
                control = list(adapt_delta = 0.9999),
                chains = 4, cores = 4, iter = 2000)

median.pars.cass1 <- summary(cassfit1,pars = "beta")$summary[,6]


newdf <- data.frame(actual = rep(beta.70,5),
                    fitted = c(grplassofit1$fit$beta[,8],
                               median.pars.hs1,
                               coef(lmfit1)[-1],
                               coef(glmfit1)[-1],
                               median.pars.cass1))
newdf <- data.frame(melt(newdf),class = rep(rep(c("SGL","Horseshoe","OLS","LASSO","LN-CASS"),each = 70),2))

ggplot(data = newdf, 
       aes(x = rep(1:70,10),y = value,col = variable)) + 
  geom_point(size = 1) + 
  geom_line() + 
  facet_wrap(~ class,nrow = 1,scales = "free_y") + 
  scale_color_manual(labels=c("Actual", "Fitted"),
                     values = c("blue","red")) +
  labs(x = "Coefficient index", 
       y = "Value of coefficient",
       color = "")

#### p = 70, n = 100 ----

grplassofit2 <- cvSGL(data = list(x = cov.120,y = y.120),
                      index = rep(1:24,each = 5),
                      maxit = 10000)

hsfit2 <- stan_glm(y.120 ~ . ,
                   data = data.frame(cov.120,y.120),
                   prior = hs(df = 2),
                   family = gaussian(),
                   prior_intercept = normal(0,10),
                   adapt_delta = 0.9999)

median.pars.hs2 <- apply(as.matrix(hsfit2),2,median)[2:121]


lmfit2 <- lm(y.120 ~ . , data = data.frame(cov.120,y.120))

glmfit2 <- cv.glmnet(x = cov.120,y = y.120)

cassfit2 <- stan("simcasshered.stan",
                 data = list(N = 100,
                             ncov = 120,
                             y = as.numeric(y.120),
                             X = cov.120,
                             sigma_indic = 10,
                             mu_indic = 0,
                             ngroup = 24),
                 control = list(adapt_delta = 0.9999),
                 chains = 4, cores = 4, iter = 10000)

median.pars.cass2 <- summary(cassfit2,pars = "beta")$summary[,6]


newdf <- data.frame(actual = rep(beta.120,4),
                    fitted = c(grplassofit2$fit$beta[,5],
                               median.pars.hs2,
                               coef(glmfit2)[-1],
                               median.pars.cass2))
newdf <- data.frame(melt(newdf),class = rep(rep(c("SGL","Horseshoe","LASSO","LN-CASS"),each = 120),2))

ggplot(data = newdf, 
       aes(x = rep(1:120,8),y = value,col = variable)) + 
  geom_point(size = 1) + 
  geom_line() + 
  facet_wrap(~ class,nrow = 1,scales = "free_y") + 
  scale_color_manual(labels=c("Actual", "Fitted"),
                     values = c("blue","red")) +
  labs(x = "Coefficient index", 
       y = "Value of coefficient",
       color = "")