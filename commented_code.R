library(splines)  # used for log(\sigma^2(t))
library(fields)   # Used for distance calculations
library(MASS)     # mvrnorm - multivariate normal draws
library(emulator) # quad.form - multivariate normal draws

############## map for families


rm(list = ls())

#' @data Y - vector of some quantity
#' @param sigma2 - variance that multiplies correlation matrix
#' @param log_det - log-determinant of correlation matrix
#' @param prec - precision matrix of correlation matrix
#' @return log-density of multivariate normal

my_dmvnorm = function(Y,sigma2,log_det,prec){
  n_local = length(Y) 
  -n_local/2 * log(2 * pi) - 
    1/2 * (log_det + n_local * log(sigma2)) - 
    quad.form(prec,Y) / (2 * sigma2)
}

#' @data x - distances between wavelength and knot location
#' @param v - v is the bandwidth (sd) of normal kernel convolution
#' @return design matrix for kernel convolution

kern = function(x,v){
  dnorm(abs(x)/v)/v
}


##### MCMC settings

reps = 1e4    #Number of MCMC samples returned (post-burn/post-thinning)
burn = 15e4   #Number of MCMC iterations before saving samples returned
tune = 100    #How often to tune Metropolis-Hastings candidate variance
thin = 5      #Keep every fifth MCMC sample

##### Knot spacing for gamma(t), eta, sig^2(t), and beta(t)
dx_genus = 25
dx_spat = 25
dx_sig2 = 20
dx_beta = 25


####### Load in log-reflectances (Y), scaled environment (X), and locations
Y = readRDS("reflectances.rds")
X = readRDS("envr.rds")
loc_aster = readRDS("locs.rds")

###### Load in indices of genus and locations
genus_ind = readRDS("genus_idx.rds")
loc_ind = readRDS("loc_idx.rds")

###### number of spectra, wavelengths sample, genera, locations, 
###### and environmental predictors in dataset
N_rep = ncol(Y)
N_wave = 500
N_genus = max(genus_ind)
N_loc = max(loc_ind)
P = ncol(X)



##### ##### ##### ##### ##### ##### ##### ##### 
##### Create K_gamma using kernel convolutions
##### ##### ##### ##### ##### ##### ##### ##### 

N_knot_genus = N_wave/dx_genus + 2
wave_knots_genus = seq(450 - dx_genus/2,950 + dx_genus/2,by = dx_genus)
### for simplicity, we use mean from posterior samples from AoAS sampling scheme
v_now_genus = apply(readRDS("v_genus.rds"),2,mean)
dif_temp_genus = outer(450:949,wave_knots_genus,"-")
K_genus = kern(dif_temp_genus,v_now_genus)

##### ##### ##### ##### ##### ##### ##### ##### 
##### Create K_eta using kernel convolutions
##### ##### ##### ##### ##### ##### ##### ##### 

N_knot_spat = N_wave/dx_spat + 2
wave_knots_spat = seq(450 - dx_spat/2,950 + dx_spat/2,by = dx_spat)
v_now_spat = dx_spat ### for simplicity, we fix the bandwidth
dif_temp_spat = outer(450:949,wave_knots_spat,"-")
K_spat = kern(dif_temp_spat,v_now_spat)

##### ##### ##### ##### ##### ##### ##### ##### 
##### Create K_beta using kernel convolutions
##### ##### ##### ##### ##### ##### ##### ##### 

N_knot_beta = N_wave/dx_beta + 2
wave_knots_beta = seq(450 - dx_beta/2,950 + dx_beta/2,by = dx_beta)
v_now_beta = dx_beta ### for simplicity, we fix the bandwidth
dif_temp_beta = outer(450:949,wave_knots_beta,"-")
K_beta = kern(dif_temp_beta,v_now_beta)

##### ##### ##### ##### ##### ##### ##### ##### 
##### Create K_sigma using kernel convolutions
##### ##### ##### ##### ##### ##### ##### ##### 

wave_knots_sig2 = seq(450 - dx_sig2/2,950 + dx_sig2/2,by = dx_sig2)
K_sig2 = cbind(1,bs(450:949,knots = wave_knots_sig2,degree = 1))
n_par_sig2 = ncol(K_sig2)

## Because splines are sparse, we leverage this in sampling
## index the observations that aren't connected with each kernel 
ind_sig2_spline = as.list(apply(K_sig2,2,function(x){which(x != 0)}))

##### ##### ##### ##### ##### ##### 
##### low-rank representation of eta
##### ##### ##### ##### ##### ##### 

### number of latent GPs, chosen through cross-validation in AoAS paper
r = 10    

### fixed spatial decay parameters
phi_spat_now = c(1,1/2,1/3,1/5,1/10,1/20,1/25,1/50,1/100,1/200)

### Distances in km between sites
dd = rdist(loc_aster)

### list of correlation matrices
R = lapply(1:r,function(x){
  exp(- dd * phi_spat_now[r])
})

### Invert correlation matrices
R_chol = lapply(1:r,function(x){t(chol(R[[x]]))})
log_det_R = sapply(R_chol,function(x){2*sum(log(diag(x)))})
inv_R = lapply(1:r,function(x){chol2inv(t(R_chol[[x]]))})


### Calculate the conditional correlation - proj is intermediate step.

proj = lapply(1:N_loc,function(x){
  sapply(1:r, function(qq){
    R[[qq]][x,-x] %*% solve(R[[qq]][-x,-x])
  })
})

cond_cor = lapply(1:N_loc,function(x){
  sapply(1:r, function(qq){
    c(R[[qq]][x,x] - proj[[x]][,qq] %*% R[[qq]][-x,x])
  })
})




### These are not necessary if bandwidths
###  aren't treated as unknowns

# genus_chol_var = t(chol(exp(-rdist(wave_knots_genus) / 50 )))
# log_det_genus = 2 * sum(log(diag(genus_chol_var)))
# inv_R_genus = chol2inv(t(genus_chol_var))
# sum_R_genus = sum(inv_R_genus)
# R_genus_colsum = apply(inv_R_genus,2,sum)
# 
# spat_chol_var = t(chol(exp(- rdist(wave_knots_spat) / 50 )))
# log_det_spat = 2 * sum(log(diag(spat_chol_var)))
# inv_R_spat = chol2inv(t(spat_chol_var))
# sum_R_spat = sum(inv_R_spat)
# R_spat_colsum = apply(inv_R_spat,2,sum)


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#####. Create indices needed for Gibbs sampling
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

#### number observations
obs_ind = 1:N_rep

#### List of observation indices at each site
loc_idx_list = lapply(1:N_loc,function(x){ which(loc_ind == x)})

#### List of locations indices for each genus
loc_by_genus_unique = lapply(1:N_genus,function(x){ unique(loc_ind[genus_ind == x]) })

#### List of locations counting from 1, including duplication
loc_by_genus_reind = lapply(1:N_genus,function(x){ as.numeric(as.factor(loc_ind[genus_ind == x])) })
N_genus_loc_count = sapply(loc_by_genus_reind,table) ### associated counts

#### List of observations indices for each genus
obs_by_genus = lapply(1:N_genus,function(x){ obs_ind[genus_ind == x] })

#### Number of unique genus/locations combos for alpha_i(s)
N_genus_loc = sum(sapply(loc_by_genus_unique,length))

#### indices used to map alpha_i(s) to observations
genus_ind_start = c(0,cumsum(sapply(loc_by_genus_reind,max)[-N_genus]))
alp_is_ind_list =  lapply(1:N_genus,function(x){ unique(genus_ind_start[x] + loc_by_genus_reind[[x]])} )
alp_is_ind_vec =  unlist(lapply(1:N_genus,function(x){ genus_ind_start[x] + loc_by_genus_reind[[x]] }))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Memory allocated for samples - all parameter names match those in 
### GMP and White et al. (2021).
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### genus-varying spatial intercept
alpha_si_now =  rep(0,max(alp_is_ind_vec)); alpha_si = matrix(0,reps,max(alp_is_ind_vec))


sig2_alp_s_now = 1; sig2_alp_s = numeric(reps)
phi_alp_s_now = 1/50; phi_alp_s = numeric(reps)


### ### ### ### ### ### ### ### 
### Create all genus-specific spatial random effects 
### ### ### ### ### ### ### ### 
### 1/50 is starting value

temp = exp(- dd * phi_alp_s_now)  ### correlation matrix with decay parameter 

Spat_Genus_R= lapply(1:N_genus,function(x){
  temp[loc_by_genus_unique[[x]],loc_by_genus_unique[[x]]]
}) 

Spat_Genus_chol = lapply(1:N_genus,function(x){
  chol(Spat_Genus_R[[x]])
}) 

Spat_Genus_inv_R = lapply(1:N_genus,function(x){
  chol2inv(Spat_Genus_chol[[x]])
}) 

Spat_Genus_log_det_R = lapply(1:N_genus,function(x){
  2*sum(log(diag(Spat_Genus_chol[[x]])))
}) 


delta_i_now = tapply(apply(Y,2,mean),genus_ind,mean); delta_i = matrix(0,reps,N_genus) 
mu_delta_now = mean(delta_i_now); mu_delta = numeric(reps) 

beta_now = matrix(0,P,N_knot_beta); beta = array(0,c(reps,P,N_knot_beta))

gamma_now = numeric(N_knot_genus); gamma = array(0,c(reps,N_knot_genus))

z_now = matrix(0,r,N_loc); z = array(0,c(reps,r,N_loc))
V_now = rep(1,r); V = array(0,c(reps,r))

A_now = matrix(1,N_knot_spat,r); A = array(0,c(reps,N_knot_spat,r))

#######

# v_genus = matrix(v_now_genus,reps,N_knot_genus)
# v_spat = numeric(reps)
# v_beta = numeric(reps)
# m_v_genus_now = mean(v_now_genus); m_v_genus = rep(m_v_genus_now,reps)
# v_v_genus_now = 1; v_v_genus = rep(1,reps)


beta_sig2_now = c(log(0.01),rep(0,n_par_sig2-1)); beta_sig2 = matrix(0,reps,n_par_sig2)
sig2_now = c(exp(K_sig2 %*% beta_sig2_now))

sig2_A_now = 1; sig2_A = numeric(reps)
sig2_bet_now = 1; sig2_bet = numeric(reps)
sig2_del_now = 1; sig2_del = numeric(reps)
sig2_gamma_now = 1; sig2_gamma = numeric(reps)

like_save = numeric(reps)

n_i_genus = c(table(genus_ind))
n_i_loc = c(table(loc_ind))

XtX = crossprod(X)

xb = K_beta %*% t(X %*%  beta_now)

Kgam_now = K_genus %*% gamma_now

KA = (K_spat %*% A_now)
KZ_now = KA %*% z_now


KtK_spat = quad.form(Matrix::diag(1/sig2_now),K_spat)
KAtKA_spat = quad.form(Matrix::diag(1/sig2_now),KA)
KtK_genus= quad.form(Matrix::diag(1/sig2_now),K_genus)
KtK_beta = quad.form(Matrix::diag(1/sig2_now),K_beta)

S_inv_beta = diag(P)/10

# cand_v_genus = rep(0.05,N_knot_genus)
# count_v_genus = rep(0,N_knot_genus)
# count_v_spat = 0
# cand_v_spat = 1e-4
# count_v_beta = 0
# cand_v_beta = 1e-2

cand_v_sig2 = rep(0.2,n_par_sig2)
count_v_sig2 = rep(0,n_par_sig2)

cand_phi_tau_alp = 0.1
count_phi_alp = 0

# e.g., use 'sweep' to express a given matrix in terms of distance from 
# the respective column mean
st = proc.time()

for(i in 2:(reps*thin + burn)){
  
  sig2_inv_sum = sum(1/sig2_now)
  
  ######### delta_i -- Gibbs update using posterior conditional distribution
  
  mu_adj = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)] + xb,2,alpha_si_now[alp_is_ind_vec],"+")
  
  
  m_del_i = tapply(apply(sweep(Y - mu_adj,1,sig2_now,"/"),2,sum),genus_ind,sum) + 
    mu_delta_now/sig2_del_now
  
  v_del_i = 1 / ( n_i_genus* sum(1/ sig2_now) + 1/sig2_del_now)
  
  delta_i_now = rnorm(N_genus,m_del_i*v_del_i,sqrt(v_del_i))
  
  ######### alpha_i(s) -- Gibbs update using posterior conditional distribution
  
  mu_adj = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)]+ xb ,2,delta_i_now[genus_ind],"+")
  
  
  m_si = tapply(apply(sweep(Y - mu_adj,1,sig2_now,"/"),2,sum),alp_is_ind_vec,sum)
  
  
  for(j in 1:N_genus){
    if(length(N_genus_loc_count[[j]]) == 1){
      v_si = 1/(c(Spat_Genus_inv_R[[j]])/sig2_alp_s_now  + N_genus_loc_count[[j]]* sig2_inv_sum)
      alpha_si_now[alp_is_ind_list[[j]]] = rnorm(1, v_si * m_si[alp_is_ind_list[[j]]],sqrt(v_si))
    }else{
      v_si = solve(Spat_Genus_inv_R[[j]]/sig2_alp_s_now  + diag(N_genus_loc_count[[j]])* sig2_inv_sum)
      alpha_si_now[alp_is_ind_list[[j]]] = mvrnorm(1, v_si %*% m_si[alp_is_ind_list[[j]]],v_si)
    }
    
  }
  
  ######### \sigma^2_{alpha_i(s)} -- Gibbs update using posterior conditional distribution
  
  alp_quad =sum(sapply(1:N_genus,function(j){
    quad.form(Spat_Genus_inv_R[[j]],alpha_si_now[alp_is_ind_list[[j]]])
  }))
  
  a_now = 3 + N_genus_loc/2
  b_now = 2 + alp_quad/2
  
  sig2_alp_s_now = 1/rgamma(1,a_now,b_now)
  
  
  ############## phi_\alpha -- Gibbs update using posterior conditional distribution
  
  phi_alp_s_cand = phi_alp_s_now
  phi_alp_s_cand = rlnorm(1,log(phi_alp_s_now),cand_phi_tau_alp)
  
  
  if(phi_alp_s_cand < 1 & phi_alp_s_cand > 1/100){
    
    temp_cand = exp(- dd * phi_alp_s_cand) 
    
    
    Spat_Genus_R_cand = lapply(1:N_genus,function(x){
      temp_cand[loc_by_genus_unique[[x]],loc_by_genus_unique[[x]]]
    }) 
    
    Spat_Genus_chol_cand = lapply(1:N_genus,function(x){
      chol(Spat_Genus_R_cand[[x]])
    }) 
    
    Spat_Genus_inv_R_cand = lapply(1:N_genus,function(x){
      chol2inv(Spat_Genus_chol_cand[[x]])
    }) 
    
    Spat_Genus_log_det_R_cand = lapply(1:N_genus,function(x){
      2*sum(log(diag(Spat_Genus_chol_cand[[x]])))
    }) 
    
    
    
    MH_dif = dlnorm(phi_alp_s_now,log(phi_alp_s_cand),cand_phi_tau_alp,log = TRUE) -
      dlnorm(phi_alp_s_cand,log(phi_alp_s_now),cand_phi_tau_alp,log = TRUE)
    
    prior_dif = sum(sapply(1:N_genus,function(x){ 
      my_dmvnorm(alpha_si_now[alp_is_ind_list[[x]]],sig2_alp_s_now,
                 Spat_Genus_log_det_R_cand[[x]],Spat_Genus_inv_R_cand[[x]])})) -
      sum(sapply(1:N_genus,function(x){ 
        my_dmvnorm(alpha_si_now[alp_is_ind_list[[x]]],sig2_alp_s_now,
                   Spat_Genus_log_det_R[[x]],Spat_Genus_inv_R[[x]])}))
    
    if(prior_dif + MH_dif > log(runif(1))){
      
      phi_alp_s_now = phi_alp_s_cand
      
      Spat_Genus_R = Spat_Genus_R_cand
      Spat_Genus_chol = Spat_Genus_chol_cand
      Spat_Genus_log_det_R = Spat_Genus_log_det_R_cand
      Spat_Genus_inv_R = Spat_Genus_inv_R_cand
      
      count_phi_alp = count_phi_alp + 1
    }
    
  }
  
  ######### mu_delta -- Gibbs update using posterior conditional distribution
  
  m_del = 0 / 100 + sum(delta_i_now) / sig2_del_now
  
  v_del = 1/(1 / 100 + N_genus / sig2_del_now)
  
  mu_delta_now = rnorm(1,m_del*v_del,sqrt(v_del))
  
  
  ######### sig2_delta -- Gibbs update using posterior conditional distribution
  
  a_sig_del = 3 + N_genus/2
  b_sig_del = 0.2 + sum((delta_i_now - mu_delta_now)^2)/2
  
  sig2_del_now = 1/rgamma(1,a_sig_del,b_sig_del)
  
  
  ######## beta -- Gibbs update using posterior conditional distribution
  
  mu_adj = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)],2,
                 delta_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec] ,"+")
  
  mu_temp = sweep(Y - mu_adj,1,sig2_now,"/")
  
  m_bet  = c(t(K_beta) %*% mu_temp %*% X)
  v_bet  = solve(diag(P*N_knot_beta)/ sig2_bet_now + XtX %x% KtK_beta)
  beta_now = matrix(mvrnorm(1,v_bet %*% m_bet,v_bet),P,N_knot_beta,byrow = TRUE)
  xb = K_beta %*% t(X %*%  beta_now)
  
  ######## \sigma^2_\beta -- Gibbs update using posterior conditional distribution
  
  a_bet = 3 + N_knot_beta * P / 2
  b_bet = 2 + 0.5 * sum(beta_now^2)
  sig2_bet_now = 1/rgamma(1,a_bet,b_bet)
  
  ########   ########   ########   ########   ########   ######## 
  ######## update gamma, gamma_i, Sig2_gamma    
  ########   ########   ########   ########   ########   ######## 
  
  ######## gamma -- Gibbs update using posterior conditional distribution
  
  mu_adj = sweep(KZ_now[,loc_ind]+ xb,2, delta_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec] ,"+")
  
  
  m_gam = t(K_genus) %*% apply(sweep(Y - mu_adj,1,sig2_now,"/"),1,sum)
  v_gam = solve(diag(N_knot_genus)/ sig2_gamma_now + N_rep * KtK_genus)
  gamma_now = mvrnorm(1,v_gam %*% m_gam ,v_gam)
  
  Kgam_now = K_genus %*% gamma_now
  
  ######## \sigma^2_\gamma -- Gibbs update using posterior conditional distribution
  
  a_gam = 3 + N_knot_genus/2
  b_gam = 2 + 0.5 * sum(gamma_now^2)
  sig2_gamma_now = 1/rgamma(1,a_gam,b_gam)
  
  
  ########   ########   ########   ########   ########   ######## 
  ######## update Z and A that make up eta(s,t)
  ########   ########   ########   ########   ########   ######## 
  
  ##### Update Z -- Gibbs update using posterior conditional distribution
  
  
  mu_adj = sweep(Kgam_now[,rep(1,N_rep)] + xb,2,delta_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec],"+")
  mu_temp = sweep(Y - mu_adj,1,sig2_now,"/")
  
  for(j in 1:N_loc){
    
    v_z = solve( diag(1/ cond_cor[[j]]) + n_i_loc[j] * KAtKA_spat)
    
    if(n_i_loc[j] == 1){
      m_z = c(mu_temp[,loc_idx_list[[j]]] %*% KA) + 
        apply(proj[[j]] * t(z_now[,-j]),2,sum) / cond_cor[[j]]
    } else{
      m_z = c(t(apply(mu_temp[,loc_idx_list[[j]]], 1,sum) %*% KA)) + 
        apply(proj[[j]] * t(z_now[,-j]),2,sum) / cond_cor[[j]]
    }
    
    
    z_now[,j] = c(mvrnorm(1,v_z %*% m_z,v_z))
    
  }
  
  ##### Update A -- Gibbs update using posterior conditional distribution
  
  for(j in 1:r){
    
    if(r == 2){
      
      mu_adj = sweep(Kgam_now[,rep(1,N_rep)] + xb + outer(KA[,-j], z_now[-j,loc_ind]),2,
                     delta_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec],"+")
      mu_temp = sweep(Y - mu_adj,1,sig2_now,"/")
      
      v_r = solve( diag(N_knot_spat) / sig2_A_now + KtK_spat * sum(z_now[j,loc_ind]^2) )
      m_r = t(K_spat) %*% apply(sweep(mu_temp,2,z_now[j,loc_ind],"*"),1,sum)
      
    } else{
      
      mu_adj = sweep(Kgam_now[,rep(1,N_rep)]  + xb + KA[,-j] %*% z_now[-j,loc_ind],2,
                     delta_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec],"+")
      mu_temp = sweep(Y - mu_adj,1,sig2_now,"/")
      
      v_r = solve( diag(N_knot_spat) / sig2_A_now + KtK_spat * sum(z_now[j,loc_ind]^2) )
      m_r = t(K_spat) %*% apply(sweep(mu_temp,2,z_now[j,loc_ind],"*"),1,sum)
      
      
    }
    
    A_now[,j] = c(mvrnorm(1,v_r %*% m_r,v_r))
    
    KA = (K_spat %*% A_now)
    
  }
  
  KZ_now = KA %*% z_now
  
  mu_now = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)]+ xb,2,
                 delta_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec] ,"+")
  
  ##### Update \sigma^2_A -- Gibbs update using posterior conditional distribution
  
  a_A = 11 + N_knot_spat * r / 2
  b_A = 10 + sum(A_now^2)/2
  sig2_A_now = 1/rgamma(1,a_A,b_A)
  
  ########   ########   ########   ########   ########   ######## 
  ######## update \sigma^2(t) using metropolis hastings
  ########   ########   ########   ########   ########   ######## 
  
  for(j in 1:n_par_sig2){
    
    beta_sig2_cand = beta_sig2_now
    sig2_cand = sig2_now
    
    beta_sig2_cand[j] = rnorm(1,beta_sig2_now[j],cand_v_sig2[j])
    
    sig2_cand = c(exp(K_sig2 %*% beta_sig2_cand))
    
    if(j == 1){
      prior_dif = dnorm(beta_sig2_cand[j],0,100,log = TRUE) - 
        dnorm(beta_sig2_now[j],0,100,log = TRUE)
    } else{
      prior_dif = dnorm(beta_sig2_cand[j],0,3,log = TRUE) - 
        dnorm(beta_sig2_now[j],0,3,log = TRUE)
    }
    
    like_dif = sum(dnorm(c(Y[ind_sig2_spline[[j]],]),c(mu_now[ind_sig2_spline[[j]],]),
                         rep(sqrt(sig2_cand[ind_sig2_spline[[j]]]),times = N_rep),log = TRUE)) - 
      sum(dnorm(c(Y[ind_sig2_spline[[j]],]),c(mu_now[ind_sig2_spline[[j]],]),
                rep(sqrt(sig2_now[ind_sig2_spline[[j]]]),times = N_rep),log = TRUE)) 
    
    if(like_dif + prior_dif > log(runif(1))){
      
      sig2_now = sig2_cand
      beta_sig2_now = beta_sig2_cand
      
      count_v_sig2[j] = count_v_sig2[j] + 1
    }
    
    
  }
  
  mod_var = rep(sqrt(sig2_now),times = N_rep)
  
  
  KtK_spat = quad.form(Matrix::diag(1/sig2_now),K_spat)
  KAtKA_spat = quad.form(Matrix::diag(1/sig2_now),KA)
  KtK_genus= quad.form(Matrix::diag(1/sig2_now),K_genus)
  KtK_beta = quad.form(Matrix::diag(1/sig2_now),K_beta)
  
  
  if(i  %% tune == 0){
    
    if(i < burn){
      
      ## Because we aren't updating bandwidths here, several of these are 
      # acc_v_genus = count_v_genus / tune; count_v_genus = rep(0,N_knot_genus)
      # acc_v_spat = count_v_spat / tune; count_v_spat = rep(0,1)
      # acc_v_beta = count_v_beta / tune; count_v_beta = rep(0,1)
      
      acc_v_sig2 = count_v_sig2 / tune; count_v_sig2 = rep(0,n_par_sig2)
      acc_phi_alp = count_phi_alp / tune; count_phi_alp = 0
      
      # cand_v_genus = ifelse(acc_v_genus > 0.6, cand_v_genus*2,
      #                       ifelse(acc_v_genus < 0.2, cand_v_genus/3, cand_v_genus) )
      # cand_v_spat = ifelse(acc_v_spat > 0.6, cand_v_spat*2,
      #                      ifelse(acc_v_spat < 0.2, cand_v_spat/3, cand_v_spat) )
      # cand_v_beta = ifelse(acc_v_beta > 0.6, cand_v_beta*2,
      #                      ifelse(acc_v_beta < 0.2, cand_v_beta/3, cand_v_beta) )
      cand_v_sig2 = ifelse(acc_v_sig2 > 0.6, cand_v_sig2*2,
                           ifelse(acc_v_sig2 < 0.2, cand_v_sig2/3, cand_v_sig2) )
      
      
      cand_phi_tau_alp = ifelse(acc_phi_alp > 0.6, cand_phi_tau_alp*2,
                                ifelse(acc_phi_alp < 0.2, cand_phi_tau_alp/3, 
                                       cand_phi_tau_alp) )
      
    } 
    
    time_its <- (proc.time() - st)[3] / (i - 1)
    time_used <- round((proc.time() - st)[3]/(60),digits=4)
    time_left <- round(time_its * (reps*thin + burn - i )/(60),digits=4)
    cat("\r", i, " of ", reps*thin + burn,"||| Time left: ",floor(time_left/60),
        " hours",time_left%%60," minutes |||| like = ", 
        sum(dnorm(c(Y),c(mu_now),mod_var,log = TRUE))) 
    flush.console()
    
    
  }
  
  
  if(i > burn & i %% thin == 0){
    ### save parameters
    
    
    alpha_si[(i - burn)/thin,] = alpha_si_now
    sig2_alp_s[(i - burn)/thin] = sig2_alp_s_now
    
    delta_i[(i - burn)/thin,] = delta_i_now
    mu_delta[(i - burn)/thin] = mu_delta_now
    beta[(i - burn)/thin,,] = beta_now
    gamma[(i - burn)/thin,] = gamma_now
    A[(i - burn)/thin,,] = A_now
    z[(i - burn)/thin,,] = z_now
    phi_alp_s[(i - burn)/thin] = phi_alp_s_now
    
    beta_sig2[(i - burn)/thin,] = beta_sig2_now
    sig2_del[(i - burn)/thin] = sig2_del_now
    sig2_gamma[(i - burn)/thin] = sig2_gamma_now
    sig2_A[(i-burn)/thin] = sig2_A_now
    
    like_save[(i - burn)/thin] =  sum(dnorm(c(Y),c(mu_now),mod_var,log = TRUE))
  }
  
  
  
}


# save.image("aster_simplified.RData")

