## -----------------------------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
summary(lm.D9)$coef

## -----------------------------------------------------------------------------
library(ggplot2)
head(iris)

ggplot(iris, aes(x = Petal.Length, fill = Species)) +
  geom_histogram(binwidth = 0.2) +
  labs(title = "Petal Length Distribution", x = "Petal Length", y = "Frequency")


## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(114514)

## -----------------------------------------------------------------------------
my_sample <- function(x,size=length(x),prob=rep(1,length(x))) {
  # 防止输入概率和不为1 进行”归一化“
  prob <- prob/sum(prob)
  # 逆变换法
  cp <- cumsum(prob)
  U = runif(size)
  r <- x[findInterval(U,cp)+1]
  return(r)
}

## -----------------------------------------------------------------------------
# Exp 1
table(my_sample(1:3,10000))
# Exp 2
table(my_sample(1:3,10000,prob = c(.2,.3,.5)))

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
u[u<.5] <- log(2*u[u<.5])
u[u>=.5] <- -log(2-2*u[u>=.5])
hist(u, prob = TRUE, main = expression(f(x)==0.5*exp(-abs(x))),xlim = c(-5,5))
y <- seq(-5, 5, .01)
lines(y, exp(-abs(y))/2)

## -----------------------------------------------------------------------------
gamma <- function(x){
  return(factorial(x-1))
}

beta_ar <- function(n, a, b){
  k = 0
  y = numeric(n) 
  #Beta(a,b)的密度函数
  fbeta <- function(x){
    return(gamma(a+b)/gamma(a)/gamma(b)*x^(a-1)*(1-x)^(b-1))
  }
  x0 = as.numeric(optimize(fbeta,lower = 0, upper = 1, maximum = TRUE)[[1]])
  while(k < n){
    u = runif(1)
    x = runif(1)
    if((fbeta(x)/(fbeta(x0)+0.5)) > u){
      k = k+1
      y[k] = x
    }
  }
  return(y)
}

## -----------------------------------------------------------------------------
# set.seed(114514)
n = 1000
x = beta_ar(n, 3, 2)
y = seq(0, 1, 0.01)
hist(x, probability = TRUE, main = "Beta(3,2)", ylim = c(0,2))
lines(y, dbeta(y, 3, 2))

## -----------------------------------------------------------------------------
n <- 1000
DG_simulation <- function(n){
  k <- 0
  while(k<n){
    k <- k+1
    u1 <- runif(1,-1,1)
    u2 <- runif(1,-1,1)
    u3 <- runif(1,-1,1)
    if(abs(u3)>=abs(u2) & abs(u3)>=abs(u2)){
      x[k] <- u2
    } else {
      x[k] <- u3
    }    
  }
  return(x)
  }
x <- DG_simulation(10000)
y = seq(-1, 1, 0.01)
hist(x,main = expression(f(x)==0.75*(1-x^2)),probability = TRUE)
lines(y, 0.75*(1-y^2))


## -----------------------------------------------------------------------------
rho <- c(0.5, 0.8, 1)
for(i in rho){
  d <- 1
  l <- d*i
  m <- 1e6
  pihat <- c()
  for (j in 1:100) {
    X <- runif(m,0,d/2)
    Y <- runif(m,0,pi/2)
    pihat[j] <- 2*l/d/mean(l/2*sin(Y)>X)
  }
  
  cat("When rho = ",i,"Var is ",var(pihat),"\n")
  
}

## -----------------------------------------------------------------------------
MC_anti <- function(n, anti=FALSE){
  u = runif(n/2)
  if(anti) {
    v = 1-u
    return((mean(exp(u))+mean(exp(v)))/2)
    }
  else {
    v = runif(n/2)
    u = c(u,v)
    return(mean(exp(u)))
    }
}
set.seed(0)
m = 1e4
v1 <- v2 <- numeric(1000)
for(i in 1:1000){
 v1[i] = MC_anti(m)
 v2[i] = MC_anti(m,TRUE)
}
var1 = var(v1)
var2 = var(v2)
cat("Simple var:",var1,"\t Anti var: ",var2)
cat("\nThe variance reduction is: ",(var1-var2)/var1)

## -----------------------------------------------------------------------------
g <- function(x){
  x^2*exp(-x^2/2)/sqrt(2*pi)*(x>=1)
}
x = seq(1,8,0.01)
gs <- c(expression(g(x)),expression(f1(x)),expression(f2(x)))
par(mfrow=c(1,2))
# figure of g, f1, f2
plot(x, g(x), type="l", ylab="", ylim=c(0,0.6), lwd = 2, col=1)
lines(x, dnorm(x), lwd=2, col=2)
lines(x, dgamma(x,3,2), lwd=2, col=3)
legend("topright", legend = gs, lty=1, lwd=2, inset = 0.02,col=1:3)

# figure of g/f1, g/f2
plot(x, g(x)/dnorm(x), type="l", ylab="", ylim=c(0,5), lwd = 2, col=2)
lines(x, g(x)/dgamma(x,3,2), lwd=2, col=3)
legend("topright", legend = gs[-1], lty=1, lwd=2, inset = 0.02,col=2:3)

## -----------------------------------------------------------------------------
set.seed(0)
m = 1e4
theta <- se <- numeric(2)
# using f1
x <- rnorm(m) 
fg <- g(x) / dnorm(x)
theta[1] <- mean(fg)
se[1] <- sd(fg)
# using f2
x <- rgamma(m,3,2) 
fg <- g(x) / dgamma(x,3,2)
theta[2] <- mean(fg)
se[2] <- sd(fg)
rbind(theta, se)

## -----------------------------------------------------------------------------
se^2

## -----------------------------------------------------------------------------
# Parameters for the gamma distribution
alpha <- 3  # Shape parameter
beta <- 1   # Scale parameter

# Number of samples
N <- 100000

samples <- rgamma(N, shape = alpha, rate = beta)
estimate <- (mean((samples^2 * gamma(alpha)) / beta^alpha))

# Print the estimate
cat("Monte Carlo Estimate:", estimate, "\n")



## -----------------------------------------------------------------------------
set.seed(0)
m <- 1e6
est1 <- sd1 <- 0
g <- function(x){
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
u <- runif(m) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
est1 <- mean(fg)
sd1 <- sd(fg)

## -----------------------------------------------------------------------------
M <- 10000
k <- 5 
m <- M/k # replicates per stratum
T1 <- T2 <- numeric(k)
est2 <- numeric(2)
fj <- function(x)5*exp(-x)/(1-exp(-1)) # fj(x)
g <- function(x)exp(-x)/(1+x^2)

## -----------------------------------------------------------------------------
set.seed(0)
for(j in 1:k){
  u = runif(m)
  x = -log(1-(1-exp(-1))*(u+j-1)/5) # inverse transform method
  T1[j] = mean(g(x)/fj(x))
  T2[j] = var(g(x)/fj(x))
}
est2[1] = sum(T1) 
est2[2] = sum(T2)
round(c(est1,est2[1]),4)
round(c(sd1,sqrt(est2[2])),5)

## -----------------------------------------------------------------------------
cl <- 0.95
s <- 20
n <- 10000
true_mean <- 2  

successful_intervals <- 0

for (i in 1:n) {
  x <- rchisq(s, df = 2)
  
  sample_mean <- mean(x)
  sample_sd <- sd(x)
  
  moe <- qt((1 + cl) / 2, df = s - 1) * (sample_sd / sqrt(s))
  interval <- c(sample_mean - moe, sample_mean + moe)
  
  if (true_mean >= interval[1] && true_mean <= interval[2]) {
    successful_intervals <- successful_intervals + 1
  }
}

c_p <- successful_intervals / n

cat("Estimated Coverage Probability:", c_p, "\n")




## -----------------------------------------------------------------------------
n <- 10000
alpha <- 0.05

true_mean_chi_sq <- 1
true_mean_uniform <- 1
true_mean_exponential <- 1

rej_chi_sq <- 0
rej_uniform <- 0
rej_exponential <- 0

# Perform Monte Carlo simulations for each case
for (i in 1:n) {
  
  sample_chi_sq <- rchisq(30, df = 1)  
  sample_uniform <- runif(30, min = 0, max = 2) 
  sample_exponential <- rexp(30, rate = 1)  
  
  t_test_chi_sq <- t.test(sample_chi_sq, mu = true_mean_chi_sq, alternative = "two.sided")
  t_test_uniform <- t.test(sample_uniform, mu = true_mean_uniform, alternative = "two.sided")
  t_test_exponential <- t.test(sample_exponential, mu = true_mean_exponential, alternative = "two.sided")
  
  if (t_test_chi_sq$p.value < alpha) {
    rej_chi_sq <- rej_chi_sq + 1
  }
  
  if (t_test_uniform$p.value < alpha) {
    rej_uniform <- rej_uniform + 1
  }
  
  if (t_test_exponential$p.value < alpha) {
    rej_exponential <- rej_exponential + 1
  }
}

t1e_chi_sq <- rej_chi_sq / n
t1e_uniform <- rej_uniform / n
t1e_exponential <- rej_exponential / n

# Print the results
cat("Type I Error Rate for χ²(1) Data:", t1e_chi_sq, "\n")
cat("Type I Error Rate for Uniform(0,2) Data:", t1e_uniform, "\n")
cat("Type I Error Rate for Exponential(1) Data:", t1e_exponential, "\n")


## -----------------------------------------------------------------------------
set.seed(114514)
# bonferroni
m <- 1000
M <- 1000
alpha <- 0.1
FWER <- numeric(M)
FDR <- numeric(M)
TPR <- numeric(M)

for (i in 1:M){
  p1 <- runif(950)
  p2 <- rbeta(50,0.1,1)
  p <- c(p1,p2)
  p <- p.adjust(p,method = 'bonferroni')
  FWER[i] <- sum(p[1:950]<alpha)>0
  FDR[i] <- sum(p[1:950]<alpha)/sum(p<alpha)
  TPR[i] <- sum(p[1:950]>alpha)/(sum(p[1:950]>alpha)+sum(p[951:1000]<alpha))
}
cat("Mean FWER:", mean(FWER), "\n")
cat("Mean FDR:", mean(FDR), "\n")
cat("Mean TPR:", mean(TPR), "\n")

## -----------------------------------------------------------------------------
set.seed(114514)
# BH
m <- 1000
M <- 1000
alpha <- 0.1
FWER <- numeric(M)
FDR <- numeric(M)
TPR <- numeric(M)

for (i in 1:M){
  p1 <- runif(950)
  p2 <- rbeta(50,0.1,1)
  p <- c(p1,p2)
  p <- p.adjust(p,method = 'BH')
  FWER[i] <- sum(p[1:950]<alpha)>0
  FDR[i] <- sum(p[1:950]<alpha)/sum(p<alpha)
  TPR[i] <- sum(p[1:950]>alpha)/(sum(p[1:950]>alpha)+sum(p[951:1000]<alpha))
}
cat("Mean FWER:", mean(FWER), "\n")
cat("Mean FDR:", mean(FDR), "\n")
cat("Mean TPR:", mean(TPR), "\n")

## -----------------------------------------------------------------------------
#set.seed(114514)
lambda <- 2
n <- c(5,10,20)
B <- 1000
m <- 1000

for (j in n){
  x <- rexp(j,lambda)
  xstar <- numeric(m)
  for (i in 1:m){
    xstar[i] <- 1/mean(sample(x,B,replace = T))
  }
  cat("n=",j,"\nBootstrap method vs theoretical on Bias:",mean(xstar)-1/mean(x),lambda*j/(j-1)-1/mean(x),"\n")
  cat("Bootstrap method vs theoretical on SE:",sd(xstar),lambda*j/((j-1)*sqrt(j-2)),"\n")
}



## -----------------------------------------------------------------------------
library(bootstrap)

#set up the bootstrap
B <- 200 #number of replicates
n <- nrow(law) #sample size
R <- numeric(B) #storage for replicates
#bootstrap estimate of standard error of R
r <- cor(law$LSAT,law$GPA)
for (b in 1:B) {
#randomly select the indices
i <- sample(1:n, size = n, replace = TRUE)
LSAT <- law$LSAT[i] #i is a vector of indices
GPA <- law$GPA[i]
R[b] <- cor(LSAT, GPA)
}
#output
t <- (R-r)/sd(R)
cat("Bootstrap t Confidence Interval for Correlation:",r+quantile(t, c(0.025, 0.975))*sd(R))

## -----------------------------------------------------------------------------
data("scor",package="bootstrap")
data1 = scor
n = nrow(data1)
sigma.hat = matrix(0,5,5)
lambda.hat = numeric(5)
theta.jack = numeric(n)
for(i in 1:n){
  data_jack = data1[-i,]  # leave-one-out
  sigma.hat = (n-2)*cov(data_jack)/(n-1) # MLE of Sigma
  lambda.hat = eigen(sigma.hat)$values
  theta.jack[i] = lambda.hat[1]/sum(lambda.hat)
}

## -----------------------------------------------------------------------------
sigma = (n-1)*cov(data1)/n
lambda.hat = eigen(sigma)$values
theta.hat = lambda.hat[1]/sum(lambda.hat)
bias.jack = (n-1)*(mean(theta.jack)-theta.hat)
se.jack = sqrt((n-1)*mean((theta.jack-mean(theta.hat))^2))
round(c(bias.jack=bias.jack,se.jack=se.jack),3)

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic) 
e1 <- e2 <- e3 <- e4 <- matrix(0,n,n) # store the squared prediction errors

# fit models on leave-two-out samples
for (i in 1:(n-1)){
  for (j in (i+1):n){
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2]*chemical[c(i,j)]
    e1[i,j] <- mean((magnetic[c(i,j)] - yhat1)^2)
  
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1]+J2$coef[2]*chemical[c(i,j)]+J2$coef[3]*chemical[c(i,j)]^2
    e2[i,j] <- mean((magnetic[c(i,j)] - yhat2)^2)
  
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(i,j)]
    yhat3 <- exp(logyhat3)
    e3[i,j] <- mean((magnetic[c(i,j)] - yhat3)^2)
    
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(i,j)])
    yhat4 <- exp(logyhat4)
    e4[i,j] <- mean((magnetic[c(i,j)] - yhat4)^2)
  }
}

## -----------------------------------------------------------------------------
2*n/(n-1)*c(mean(e1), mean(e2), mean(e3), mean(e4))

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}

n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)

R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:50
n <- length(x)
m <- length(y)
reps <- numeric(R) #storage for replicates
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = n, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
reps[i] <- count5test(x1,y1)
}
mean(reps)

## -----------------------------------------------------------------------------
fun <- function(N,b1,b2,b3,f0){
  x1 <- rpois(N,1)
  x2 <- rexp(N,1)
  x3 <- rbinom(N,1,0.5)
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3); p <- 1/(1+tmp)
    mean(p) - f0
  }
  solution <- uniroot(g,c(-20,0))
  solution$root
}
f0 <- c(0.1,0.01,0.001,0.0001)
a <- c()
for (i in 1:4){
  a[i] <- fun(1e6,0,1,-1,f0[i])
}
plot(-log(f0),a,type = "b")

## -----------------------------------------------------------------------------
f <- function(x){
  return(exp(-abs(x))/2)
}
invF <- function(x){
  return(-log(2*(1-x)))
}

rw.Laplace <- function(sigma, x0, N){
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if (u[i] <= (f(y)/f(x[i-1])))
    x[i] = y else {
    x[i] = x[i-1]
    k = k + 1
  }
}
return(list(x=x, k=k))
}

Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
  psi = as.matrix(psi)
  n = ncol(psi)
  k = nrow(psi)
  psi.means = rowMeans(psi) #row means
  B = n * var(psi.means) #between variance est.
  psi.w = apply(psi, 1, "var") #within variances
  W = mean(psi.w) #within est.
  v.hat = W*(n-1)/n + (B/n) #upper variance est.
  r.hat = v.hat / W #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
set.seed(1234)
n = 12000 #length of chains
k = 4 #number of chains to generate
x0 = c(-10,-5,5,10) #different initial values
sigma = c(0.5, 2, 5) #different variance
b = 1000 #burn-in length
X1 <- X2 <- X3 <- matrix(0,k,n)
rej = matrix(0,3,4)
refline = c(-invF(0.975), invF(0.975))

for(i in 1:k){
  rw1 = rw.Laplace(sigma[1], x0[i], n)
  rw2 = rw.Laplace(sigma[2], x0[i], n)
  rw3 = rw.Laplace(sigma[3], x0[i], n)
  X1[i,]=rw1$x;X2[i,]=rw2$x;X3[i,]=rw3$x
  rej[,i] = c(rw1$k,rw2$k,rw3$k)
}

# compute diagnostic statistics
psi1 <- t(apply(X1, 1, cumsum))
psi2 <- t(apply(X2, 1, cumsum))
psi3 <- t(apply(X3, 1, cumsum))

for (i in 1:k){
  psi1[i,] = psi1[i,] / (1:ncol(psi1))
  psi2[i,] = psi2[i,] / (1:ncol(psi2))
  psi3[i,] = psi3[i,] / (1:ncol(psi3))
}
print(c(Gelman.Rubin(psi1),Gelman.Rubin(psi2),Gelman.Rubin(psi3)))

## -----------------------------------------------------------------------------
1-apply(rej, 1, mean)/n

## -----------------------------------------------------------------------------
Gibbs.binorm <- function(N,mu1=0,mu2=0,sigma1=1,sigma2=1,rho=0.9,x0=c(mu1,mu2)){
  X = matrix(0, N, 2) #the chain, a bivariate sample
  s1 = sqrt(1-rho^2)*sigma1
  s2 = sqrt(1-rho^2)*sigma2
  X[1,] = x0
  for (i in 2:N) {
    y1 = X[i-1, 2]
    m1 = mu1 + rho * (y1 - mu2) * sigma1/sigma2
    X[i, 1] = rnorm(1, m1, s1)
    x1 = X[i, 1]
    m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] = rnorm(1, m2, s2)
  }
  return(X)
}

## -----------------------------------------------------------------------------
set.seed(1234)
k = 4 #number of chains to generate
N = 10000 #length of chain
X <- matrix(0,N,2*k)
Xt <- Yt <- matrix(0,k,N)
x0 = matrix(rep(c(-1,-0.5,0.5,1),2),4) #initial values
psix = psiy = matrix(0,k,N)

for(i in 1:k){
  X[,(2*i-1):(2*i)] = Gibbs.binorm(N,x0=x0[i,])
  Xt[i,] = t(X[,2*i-1])
  Yt[i,] = t(X[,2*i])
}
psix = t(apply(Xt, 1, cumsum))
psiy = t(apply(Yt, 1, cumsum))
for (i in 1:k){
  psix[i,] = psix[i,] / (1:ncol(psix))
  psiy[i,] = psiy[i,] / (1:ncol(psiy))
}
print(c(Gelman.Rubin(psix),Gelman.Rubin(psiy)))

## ----out.width='50%', out.height='50%'----------------------------------------
b = 1000 #burn-in length
rhatx = rhaty = rep(0,N)
for (j in (b+1):N){
  rhatx[j] <- Gelman.Rubin(psix[,1:j])
  rhaty[j] <- Gelman.Rubin(psiy[,1:j])
}


## -----------------------------------------------------------------------------
x <- X[(b+1):N, 1:2]
plot(x, main="", cex=.5, xlab=bquote(X[t]),
ylab=bquote(Y[t]), ylim=range(x[,2]))

## -----------------------------------------------------------------------------
fit = lm(x[,2] ~ x[,1]) #linear model Y=a+bX
summary(fit)

## ----out.width='50%', out.height='50%'----------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
f <- function(x, sigma) {
if (any(x < 0)) return (0)
stopifnot(sigma > 0)
return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

m <- 10000
sigma <- 4
x <- numeric(m)
x[1] <- rchisq(1, df=1)
k <- 0
u <- runif(m)
for (i in 2:m) {
xt <- x[i-1]
y <- rchisq(1, df = xt)
num <- f(y, sigma) * dchisq(xt, df = y)
den <- f(xt, sigma) * dchisq(y, df = xt)
if (u[i] <= num/den) x[i] <- y else {
x[i] <- xt
k <- k+1 #y is rejected
}
}

index <- 5000:5500
y1 <- x[index]
plot(index, y1, type="l", main="", ylab="x")

## -----------------------------------------------------------------------------
x0 <- x
x <- numeric(m)
x[1] <- rchisq(1, df=1)
k <- 0
u <- runif(m)
for (i in 2:m) {
xt <- x[i-1]
y <- rchisq(1, df = xt)
num <- f(y, sigma) * dchisq(xt, df = y)
den <- f(xt, sigma) * dchisq(y, df = xt)
if (u[i] <= num/den) x[i] <- y else {
x[i] <- xt
k <- k+1 #y is rejected
}
}

index <- 5000:5500
y1 <- x[index]
plot(index, y1, type="l", main="", ylab="x")

## -----------------------------------------------------------------------------
library(coda)
summary(as.mcmc(x))

## -----------------------------------------------------------------------------
mcmc_obj <- mcmc.list(mcmc(x0),mcmc(x))

# 查看 summary
summary(mcmc_obj)
gelman.diag(mcmc_obj)
gelman.plot(mcmc_obj)

## -----------------------------------------------------------------------------
data = c(11,12,8,9,27,28,13,14,16,17,0,1,23,24,10,11,24,25,2,3)
X = matrix(data, nrow=10, byrow=TRUE)
u = X[,1]; v = X[,2]
logL <- function(lambda){
  s = sum(log(exp(-lambda*u)-exp(-lambda*v)))
  return(s)
}
res = optimize(logL, lower=0, upper=5, maximum=TRUE)
lambdahat = res[[1]]

## -----------------------------------------------------------------------------
n = 10
delta = 1
lambda = numeric(100) 
lambda[1] = 1 #initialize 
i = 1
# E-step
Econd <- function(lambda,u,v){ 
  (u*exp(-lambda*u)-v*exp(-lambda*v))/(exp(-lambda*u)-exp(-lambda*v))+1/lambda
}
while(delta > 1e-5){
  lambda[i+1] =  n/sum(Econd(lambda[i],u,v)) # M-step
  delta = abs(lambda[i+1]-lambda[i])
  i = i+1
}
lambdahat.EM = lambda[i]
round(c(lambdahat,lambdahat.EM),5)

## -----------------------------------------------------------------------------
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
library(boot) #needed for simplex function
B <- A+2

## -----------------------------------------------------------------------------
solve.game <- function(A) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
min.A <- min(A)
A <- A - min.A #so that v >= 0
max.A <- max(A)
A <- A / max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
#the ’solution’ is [x1,x2,...,xm | value of game]
#
#minimize v subject to ...
#let y strategies 1:n, with v as extra variable
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A * max.A + min.A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max.A + min.A)
soln
}

s1 <- solve.game(A)
s2 <- solve.game(B)
round(cbind(s1$x, s2$x), 7)
round(cbind(s1$v, s2$v), 7)

## -----------------------------------------------------------------------------
x = data.frame(c())
nrow(x);ncol(x)

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}
head(data.frame(lapply(iris,function(x) if(is.numeric(x)) scale01(x) else x)))

## -----------------------------------------------------------------------------
head(data.frame(lapply(iris[sapply(iris,is.numeric)],scale01)))

## -----------------------------------------------------------------------------
vapply(cars, sd, numeric(1))

## -----------------------------------------------------------------------------
vapply(iris[vapply(iris, is.numeric, logical(1))],
       sd, 
       numeric(1))

## -----------------------------------------------------------------------------
library(microbenchmark)
library(Rcpp)

sourceCpp("../src/GibbsC.cpp")

n <- 10
a <- 2
b <- 3
num_samples <- 1000

gibbs_sampler_r <- function(num_samples, n, a, b) {
  samples <- matrix(0, nrow = num_samples, ncol = 2)
  
  x <- 0
  y <- 0.5
  
  for (i in 1:num_samples) {
    x <- rbinom(1, n, y)
    y <- rbeta(1, x + a, n - x + b)
    
    samples[i, ] <- c(x, y)
  }
  
  return(samples)
}

mb <- microbenchmark(
  gibbs_sampler_r(num_samples, n, a, b),
  gibbsSampler(num_samples, n, a, b),
  times = 100
)

summary(mb)[,c(1,3,5,6)]


