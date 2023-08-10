library(Rcpp)

sourceCpp("SLCM_helper.cpp")

## parameters
K = 3 # number of skills
rho = 0 # correlation between skills
N = 1000 # number of students
J = 20 # number of questions
rep_each = 10 # number of initial values


## priors 
d_0 = rep(1,(2^K))
w_0 = c(1,100)
sigma_B = 0.5


## structure matrix (attribute profiles)
k3 = as.matrix(attr(terms(y~(x1+x2+x3)^3),"factors")[-1,]) # structure matrix
O = des_mat(rbind(rep(0),t(k3)),k3) # design matrix for each class
Q = t(k3[,c(1:3,1:3,1:3,4:6,4:6,4:5,7,7,7)]) # true Q
k3_and_0 = rbind(rep(0,3),t(k3)) # attribute profile for each class

## generate data 
gen = Data_gen(N,Q,rho,k3)

## get true beta
Beta_t = t(solve(O,qnorm(gen$theta))) # true beta 


SLCM <- function(burn_in=30000, chain_length=10000){
  
  ######## Data Generation #############
  gen = Data_gen(N,Q,rho,k3)
  # A_true = gen$a # true skill profiles 
  # Theta_true = gen$theta # true probability matrix
  Y = gen$Y
  
  
  LogL = -9999999999999999 # log likelihood
  
  for (s in 1:rep_each){
    
    ######## Initialization #############
    # Class proportions pi 
    Pi = rDirichlet(d_0)
    
    # Sparse matrix Delta 
    W_0 = rDirichlet(w_0)[1]
    w = W_0
    # Delta = matrix(1,J,(2^K-1))
    Delta = random_D(J, K, 0.1)
    
    # Coefficients matrix Beta
    Beta = matrix(0,J,2^K)
    non_zero_beta = res_norm(rep(0,sum(Delta)),rep(sigma_B,sum(Delta)),rep(1,sum(Delta)),rep(0,sum(Delta)))
    Beta[,-1][Delta>0] = non_zero_beta 
    Beta[,1] = myrnorm(rep(0,J),rep(sigma_B,J))
    
    
    ## list to save posterior samples ####
    posterior_beta = list()
    posterior_A = list()
    
    ######## Update ############
    for (t in 1:(burn_in + chain_length)){
      
      ## update piic: membership probability of each student 
      mu = O%*%t(Beta)
      Piic = update_piic(mu,Pi,Y)
      
      ## update alpha: attribute profiles
      class = update_a(Piic)
      alpha =  k3_and_0[class+1,]
      
      ## update pi: class proportions
      class_1 = c(class, 0:(2^K-1)) # in case some class has no student
      Pi = rDirichlet(as.numeric(table(class_1)) - 1 + d_0)
      
      
      ## get design matrix A 
      A = des_mat(alpha,k3)
      
      
      ## update Z
      Z = update_Z(A%*%t(Beta),Y)
      
      ## mu and sigma: mean and std of coefficients beta 
      my_par = update_mu_sigma (Z,A,Beta,sigma_B)
      Mu =  my_par$Mu
      Sigma2_p =  my_par$Sigma2_p
      
      
      ## Beta and Delta 
      D = gen_update_D(Delta,k3,Beta, Mu, sqrt(Sigma2_p), w, sigma_B, k3)
      Delta = D$Delta
      Beta = D$Beta
      
      ## update W: prior on coefficient beta's activeness 
      w = rDirichlet(w_0+c(sum(Delta),J*(2^K-1)-sum(Delta)))[1]
      
      ## save the results out
      if (t > burn_in){
        posterior_beta[[t-burn_in]] = Beta
        posterior_A[[t-burn_in]] = A
      }
    }
    
    mean_A = apply(simplify2array(posterior_A), 1:2, mean) # estimated A
    mean_beta_1= apply(simplify2array(posterior_beta), 1:2, mean) # estimated beta
    theta = pnorm(mean_A%*%t(mean_beta_1)) # estimated theta 
    
    newLogL = sum(log(theta+10^-9)*Y +log(1-theta+10^-9)*(1-Y)) # log-likelihood
    
    if(newLogL > LogL){
      LogL = newLogL
      mean_beta = mean_beta_1
    }
    
  }
  return (mean_beta)
}


result = SLCM()
