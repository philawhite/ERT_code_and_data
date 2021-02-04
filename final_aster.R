library(tidyverse)
library(splines)
library(geoR)
library(xtable)
library(fields)
library(ggplot2)
library(reshape2)
library(abind)
library(sp)
library(rgdal)
library(MASS)
library(MCMCpack)
library(emulator)
library(scoringRules)
library(parallel)

############## map for families


rm(list = ls())
dat = read.csv("SpectralReflectanceERT.csv")

########## My multivariate normal

my_dmvnorm = function(Y,sigma2,log_det,prec){
  n_local = length(Y) 
  -n_local/2 * log(2 * pi) - 
    1/2 * (log_det + n_local * log(sigma2)) - 
    quad.form(prec,Y) / (2 * sigma2)
}

dat_sub3 = dat[dat$Region %in% c("HTR","Cederberg") &
                 dat$family_GM %in% c("ASTERACEAE", "AIZOACEAE","RESTIONACEAE"),]



dim(dat_sub3)

rm(dat)

dat_sub3$genus = gsub("([A-Za-z]+).*", "\\1", dat_sub3$finalname)
dat_sub3$loc_ind = as.factor(as.numeric(factor(c(dat_sub3$latitude*dat_sub3$longitude))))


dat_sub3$loc_ind = as.factor(as.numeric(factor(c(dat_sub3$latitude*dat_sub3$longitude))))



dat_sub3$Region = factor(dat_sub3$Region)
dat_sub3$family_GM = factor(dat_sub3$family_GM)
dat_sub3$obs_ind = 1:nrow(dat_sub3)
dat_sub3$finalname = factor(dat_sub3$finalname)


table(dat_sub3$genus)

nrow(dat_sub3)


genus_use = names(table(dat_sub3$genus))[which(as.numeric(table(dat_sub3$genus)) > 0 )]


dat_sub3_genus = dat_sub3

dat_sub3_genus = dat_sub3[dat_sub3$genus %in% genus_use,] 
# 
# nrow(dat_sub3_genus)
# 
# #rm(dat_sub3)
# 
# 
# 
dat_sub3_genus$genus = factor(dat_sub3_genus$genus)
dat_sub3_genus$Region = factor(dat_sub3_genus$Region)
dat_sub3_genus$family_GM = factor(dat_sub3_genus$family_GM)
dat_sub3_genus$obs_ind = factor(1:nrow(dat_sub3_genus))
dat_sub3_genus$loc_ind = factor(dat_sub3_genus$loc_ind)
dat_sub3_genus$finalname = factor(dat_sub3_genus$finalname)


out = dat_sub3_genus %>% gather(key = "Wavelength",value = "Reflectance",
                                paste("X",450:949,sep = ""))

out = out[, c("Region","family_GM","genus","latitude","longitude","obs_ind",
              "finalname","loc_ind","Wavelength","Reflectance")]

out = out[order(out$loc_ind,out$obs_ind,out$family_GM),]
out$Wavelength = as.numeric( str_sub(out$Wavelength,2,4))
out$Wave_great_700 = as.factor(1 * (out$Wavelength > 700))

loc_unique = unique(out[,c("loc_ind","longitude","latitude")])


#######################

envr = unique(read.csv("EnvCleanERT.csv")[,-c(1,2)])

envr_loc = cbind(envr$longitude,envr$latitude)

envr_reflect = NULL 

for(i in 1:nrow(loc_unique)){
  
  idx = which(apply(envr_loc,1,function(x){all(x == loc_unique[i,-1])}))
  envr_reflect = rbind(envr_reflect,envr[idx,])
  
  out
  
}

head(envr_reflect)


out <- left_join(out, envr_reflect, by=c("longitude","latitude"))
head(out)




dat_aster = out[out$family_GM == "ASTERACEAE",]


dat_aster$genus = as.numeric(factor(dat_aster$genus))
dat_aster$Region = factor(dat_aster$Region)
dat_aster$family_GM = factor(dat_aster$family_GM)
dat_aster$obs_ind = as.numeric(factor(1:nrow(dat_aster)))
dat_aster$loc_ind = as.numeric(factor(dat_aster$loc_ind))
dat_aster$finalname = factor(dat_aster$finalname)


dim(dat_aster)[1]/500


#### For each family, we have
#### genus j with possible replication k

X_aster = unique(dat_aster[,-c(1:12)[-8]])


Y = t(matrix(log(dat_aster$Reflectance),ncol = 500,byrow = TRUE))

loc_aster = as.matrix(sp::spTransform(
  SpatialPoints(unique(cbind(dat_aster$longitude, dat_aster$latitude)), 
                proj4string=CRS("+proj=longlat")), 
  CRS("+init=epsg:32748"))@coords) / 1000



# kern = function(x,v){
#   dexp(abs(x)/v) / (2*v)
# }

kern = function(x,v){
  dnorm(abs(x)/v)/v
}

dx_genus = 25
dx_spat = 25
dx_tau = 20
dx_beta = 25

r = 10

stan_data = list(
  dx_genus = dx_genus,
  dx_spat = dx_spat,
  dx_tau = dx_tau,
  Y = Y,
  N_rep = ncol(Y),
  r = r,
  N_knot_genus = 500/dx_genus + 2,
  N_knot_spat = 500/dx_spat + 2,
  N_knot_beta = 500/dx_beta + 2,
  N_genus = length(unique(dat_aster$genus)),
  N_wave = 500,
  genus_ind = dat_aster$genus[which(1:nrow(dat_aster) %% 500 == 1)],
  loc_ind = dat_aster$loc_ind[which(1:nrow(dat_aster) %% 500 == 1)]
)



reps = 1e4
burn = 15e4
tune = 100
thin = 5


attach(stan_data)

X = scale(unique(as.matrix(dat_aster[,-c(1:12)])[,-10])[loc_ind,c(1,6,8,10)])
N_loc = max(loc_ind)
P = ncol(X)

loc_idx_list = lapply(1:N_loc,function(x){ which(loc_ind == x)})


wave_knots_genus = seq(450 - dx_genus/2,950 + dx_genus/2,by = dx_genus)
v_now_genus = dx_genus
dif_temp_genus = outer(450:949,wave_knots_genus,"-")
K_genus = kern(dif_temp_genus,v_now_genus)

wave_knots_spat = seq(450 - dx_spat/2,950 + dx_spat/2,by = dx_spat)
v_now_spat = dx_spat
dif_temp_spat = outer(450:949,wave_knots_spat,"-")
K_spat = kern(dif_temp_spat,v_now_spat)


wave_knots_beta = seq(450 - dx_beta/2,950 + dx_beta/2,by = dx_beta)
v_now_beta = dx_beta
dif_temp_beta = outer(450:949,wave_knots_beta,"-")
K_beta = kern(dif_temp_beta,v_now_beta)

wave_knots_tau = seq(450 - dx_tau/2,950 + dx_tau/2,by = dx_tau)

K_tau = cbind(1,bs(450:949,knots = wave_knots_tau,degree = 1))

ind_tau_spline = as.list(apply(K_tau,2,function(x){which(x != 0)}))
n_par_tau = ncol(K_tau)



phi_spat = matrix(1/50,reps,r); phi_spat_now = c(1,1/2,1/3,1/5,1/10,1/20,1/25,1/50,1/100,1/200)

dd = rdist(loc_aster)


R = lapply(1:r,function(x){
  exp(- dd * phi_spat_now[r])
})

R_chol = lapply(1:r,function(x){t(chol(R[[x]]))})
log_det_R = sapply(R_chol,function(x){2*sum(log(diag(x)))})
inv_R = lapply(1:r,function(x){chol2inv(t(R_chol[[x]]))})

proj = lapply(unique(dat_aster$loc_ind),function(x){
  sapply(1:r, function(qq){
    R[[qq]][x,-x] %*% solve(R[[qq]][-x,-x])
  })
})

cond_cor = lapply(unique(dat_aster$loc_ind),function(x){
  sapply(1:r, function(qq){
    c(R[[qq]][x,x] - proj[[x]][,qq] %*% R[[qq]][-x,x])
  })
})




genus_chol_var = t(chol(exp(-rdist(wave_knots_genus) / 50 )))
log_det_genus = 2 * sum(log(diag(genus_chol_var)))
inv_R_genus = chol2inv(t(genus_chol_var))
sum_R_genus = sum(inv_R_genus)
R_genus_colsum = apply(inv_R_genus,2,sum)

spat_chol_var = t(chol(exp(- rdist(wave_knots_spat) / 50 )))
log_det_spat = 2 * sum(log(diag(spat_chol_var)))
inv_R_spat = chol2inv(t(spat_chol_var))
sum_R_spat = sum(inv_R_spat)
R_spat_colsum = apply(inv_R_spat,2,sum)

# temp = lm(c(Y) ~ rep(genus_ind,each = 500) + X[rep(loc_ind,each = 500),] + 
#             K[rep(1:(N_wave),times =N_rep),])
# 
# summary(temp)

obs_ind = 1:N_rep

loc_by_genus_unique = lapply(1:N_genus,function(x){ unique(loc_ind[genus_ind == x]) })
loc_by_genus_reind = lapply(1:N_genus,function(x){ as.numeric(as.factor(loc_ind[genus_ind == x])) })
obs_by_genus = lapply(1:N_genus,function(x){ obs_ind[genus_ind == x] })

genus_ind_start = c(0,cumsum(sapply(loc_by_genus_reind,max)[-N_genus]))
alp_is_ind_list =  lapply(1:N_genus,function(x){ unique(genus_ind_start[x] + loc_by_genus_reind[[x]])} )
alp_is_ind_vec =  unlist(lapply(1:N_genus,function(x){ genus_ind_start[x] + loc_by_genus_reind[[x]] }))

N_genus_loc_count = sapply(loc_by_genus_reind,table)
N_genus_loc = sum(sapply(loc_by_genus_unique,length))

temp = exp(- dd / 50) 

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

alpha_si_now =  rep(0,max(alp_is_ind_vec)); alpha_si = matrix(0,reps,max(alp_is_ind_vec))
sig2_alp_s_now = 1; sig2_alp_s = numeric(reps)
phi_alp_s_now = 1/50; phi_alp_s = numeric(reps)

alpha_i_now = tapply(apply(Y,2,mean),genus_ind,mean); alpha_i = matrix(0,reps,N_genus) 
alpha_now = mean(alpha_i_now); alpha = numeric(reps) 

beta_now = matrix(0,P,N_knot_beta); beta = array(0,c(reps,P,N_knot_beta))

gamma_now = numeric(N_knot_genus); gamma = array(0,c(reps,N_knot_genus))

Sigma_gam_now = diag(N_knot_genus); Sigma_gam = array(0,c(reps,N_knot_genus,N_knot_genus))
Sigma_gam_inv = solve(Sigma_gam_now)

z_now = matrix(0,r,N_loc); z = array(0,c(reps,r,N_loc))
V_now = rep(1,r); V = array(0,c(reps,r))

A_now = matrix(1,N_knot_spat,r); A = array(0,c(reps,N_knot_spat,r))

#######

v_genus = matrix(v_now_genus,reps,N_knot_genus)
v_now_genus = rep(v_now_genus,N_knot_genus)

v_spat = numeric(reps)
v_beta = numeric(reps)


m_v_genus_now = mean(v_now_genus); m_v_genus = rep(m_v_genus_now,reps)
v_v_genus_now = 1; v_v_genus = rep(1,reps)


beta_tau2_now = c(log(0.01),rep(0,n_par_tau-1)); beta_tau2 = matrix(0,reps,n_par_tau)
tau2_now = c(exp(K_tau %*% beta_tau2_now))

sig2_A_now = 1; sig2_A = numeric(reps)
sig2_bet_now = 1; sig2_bet = numeric(reps)
sig2_alp_now = 1; sig2_alp = numeric(reps)
sig2_gamma_now = 1; sig2_gamma = numeric(reps)

like_save = numeric(reps)

n_i_genus = c(table(genus_ind))
n_i_loc = c(table(loc_ind))

XtX = crossprod(X)

xb = K_beta %*% t(X %*%  beta_now)

Kgam_now = K_genus %*% gamma_now

KA = (K_spat %*% A_now)
KZ_now = KA %*% z_now


KtK_spat = quad.form(Matrix::diag(1/tau2_now),K_spat)
KAtKA_spat = quad.form(Matrix::diag(1/tau2_now),KA)
KtK_genus= quad.form(Matrix::diag(1/tau2_now),K_genus)
KtK_beta = quad.form(Matrix::diag(1/tau2_now),K_beta)

S_inv_beta = diag(P)/10

cand_v_genus = rep(0.05,N_knot_genus)
count_v_genus = rep(0,N_knot_genus)
count_v_spat = 0
cand_v_spat = 1e-4

count_v_beta = 0
cand_v_beta = 1e-2

cand_v_tau = rep(0.2,n_par_tau)
count_v_tau = rep(0,n_par_tau)


cand_phi_tau = rep(0.1,r)
count_phi = rep(0,r)


cand_phi_tau_alp = 0.1
count_phi_alp = 0

# e.g., use 'sweep' to express a given matrix in terms of distance from 
# the respective column mean
st = proc.time()

for(i in 2:(reps*thin + burn)){
  
  tau2_inv_sum = sum(1/tau2_now)
  
  ######### alpha_i
  
  mu_adj = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)] + xb,2,alpha_si_now[alp_is_ind_vec],"+")
  
  
  m_alp_i = tapply(apply(sweep(Y - mu_adj,1,tau2_now,"/"),2,sum),genus_ind,sum) + 
    alpha_now/sig2_alp_now
  
  v_alp_i = 1 / ( n_i_genus* sum(1/ tau2_now) + 1/sig2_alp_now)
  
  alpha_i_now = rnorm(N_genus,m_alp_i*v_alp_i,sqrt(v_alp_i))
  
  ######### alpha_i(s)
  mu_adj = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)]+ xb ,2,alpha_i_now[genus_ind],"+")
  
  
  m_si = tapply(apply(sweep(Y - mu_adj,1,tau2_now,"/"),2,sum),alp_is_ind_vec,sum)
  
  
  for(j in 1:N_genus){
    if(length(N_genus_loc_count[[j]]) == 1){
      v_si = 1/(c(Spat_Genus_inv_R[[j]])/sig2_alp_s_now  + N_genus_loc_count[[j]]* tau2_inv_sum)
      alpha_si_now[alp_is_ind_list[[j]]] = rnorm(1, v_si * m_si[alp_is_ind_list[[j]]],sqrt(v_si))
    }else{
      v_si = solve(Spat_Genus_inv_R[[j]]/sig2_alp_s_now  + diag(N_genus_loc_count[[j]])* tau2_inv_sum)
      alpha_si_now[alp_is_ind_list[[j]]] = mvrnorm(1, v_si %*% m_si[alp_is_ind_list[[j]]],v_si)
    }
    
  }
  
  
  alp_quad =sum(sapply(1:N_genus,function(j){
    quad.form(Spat_Genus_inv_R[[j]],alpha_si_now[alp_is_ind_list[[j]]])
  }))
  
  a_now = 3 + N_genus_loc/2
  b_now = 2 + alp_quad/2
  
  sig2_alp_s_now = 1/rgamma(1,a_now,b_now)
  
  
  ############## phi_alpha
  
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

  ######### alpha
  
  m_alp = 0 / 100 + sum(alpha_i_now) / sig2_alp_now
  
  v_alp = 1/(1 / 100 + N_genus / sig2_alp_now)
  
  alpha_now = rnorm(1,m_alp*v_alp,sqrt(v_alp))
  
  
  ######### sig2_alpha
  
  a_sig_alp = 3 + N_genus/2
  b_sig_alp = 0.2 + sum((alpha_i_now - alpha_now)^2)/2
  
  sig2_alp_now = 1/rgamma(1,a_sig_alp,b_sig_alp)
  
  
  ######## beta 
  
  mu_adj = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)],2,
                 alpha_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec] ,"+")
  
  mu_temp = sweep(Y - mu_adj,1,tau2_now,"/")
  
  m_bet  = c(t(K_beta) %*% mu_temp %*% X)
  v_bet  = solve(diag(P*N_knot_beta)/ sig2_bet_now + XtX %x% KtK_beta)
  beta_now = matrix(mvrnorm(1,v_bet %*% m_bet,v_bet),P,N_knot_beta,byrow = TRUE)
  xb = K_beta %*% t(X %*%  beta_now)
  
  
  a_bet = 3 + N_knot_beta * P / 2
  b_bet = 2 + 0.5 * sum(beta_now^2)
  sig2_bet_now = 1/rgamma(1,a_bet,b_bet)
  
  ########   ########   ########   ########   ########   ######## 
  ######## updat gamma, gamma_i, Sig2_gamma    
  ########   ########   ########   ########   ########   ######## 
  
  mu_adj = sweep(KZ_now[,loc_ind]+ xb,2, alpha_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec] ,"+")
  
  
  m_gam = t(K_genus) %*% apply(sweep(Y - mu_adj,1,tau2_now,"/"),1,sum)
  v_gam = solve(diag(N_knot_genus)/ sig2_gamma_now + N_rep * KtK_genus)
  gamma_now = mvrnorm(1,v_gam %*% m_gam ,v_gam)
  
  
  Kgam_now = K_genus %*% gamma_now
  
  a_gam = 3 + N_knot_genus/2
  b_gam = 2 + 0.5 * sum(gamma_now^2)
  sig2_gamma_now = 1/rgamma(1,a_gam,b_gam)
  
  
  ########   ########   ########   ########   ########   ######## 
  ######## update Z and A
  ########   ########   ########   ########   ########   ######## 
  
  mu_adj = sweep(Kgam_now[,rep(1,N_rep)] + xb,2,alpha_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec],"+")
  mu_temp = sweep(Y - mu_adj,1,tau2_now,"/")
  
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
  
  ##### Update A
  
  
  for(j in 1:r){
    
    if(r == 2){
      
      mu_adj = sweep(Kgam_now[,rep(1,N_rep)] + xb + outer(KA[,-j], z_now[-j,loc_ind]),2,
                     alpha_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec],"+")
      mu_temp = sweep(Y - mu_adj,1,tau2_now,"/")
      
      v_r = solve( diag(N_knot_spat) / sig2_A_now + KtK_spat * sum(z_now[j,loc_ind]^2) )
      m_r = t(K_spat) %*% apply(sweep(mu_temp,2,z_now[j,loc_ind],"*"),1,sum)
      
    } else{
      
      mu_adj = sweep(Kgam_now[,rep(1,N_rep)]  + xb + KA[,-j] %*% z_now[-j,loc_ind],2,
                     alpha_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec],"+")
      mu_temp = sweep(Y - mu_adj,1,tau2_now,"/")
      
      v_r = solve( diag(N_knot_spat) / sig2_A_now + KtK_spat * sum(z_now[j,loc_ind]^2) )
      m_r = t(K_spat) %*% apply(sweep(mu_temp,2,z_now[j,loc_ind],"*"),1,sum)
      
      
    }
    
    A_now[,j] = c(mvrnorm(1,v_r %*% m_r,v_r))
    
    KA = (K_spat %*% A_now)
    
  }
  
  KZ_now = KA %*% z_now
  
  a_A = 11 + N_knot_spat * r / 2
  b_A = 10 + sum(A_now^2)/2
  sig2_A_now = 1/rgamma(1,a_A,b_A)
  
  mu_now = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)]+ xb,2,
                 alpha_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec] ,"+")
  
  ################## update K_\gamma
  
  mod_var = rep(sqrt(tau2_now),times = N_rep)
  
  for(j in 1:N_knot_genus){ #### joint sampling
    
    v_cand_genus = v_now_genus
    v_cand_genus[j] = rlnorm(1,log(v_now_genus[j]),cand_v_genus[j])
    
    if(v_cand_genus[j] < 200){
      
      K_genus_cand = K_genus
      
      MH_dif = dlnorm(v_now_genus[j],log(v_cand_genus[j]),cand_v_genus[j],log = TRUE) - 
        dlnorm(v_cand_genus[j],log(v_now_genus[j]),cand_v_genus[j],log = TRUE)
      
      
      
      K_genus_cand[,j] = kern(dif_temp_genus[,j],v_cand_genus[j])
      
      Kgam_cand = K_genus_cand %*% gamma_now
      
      mu_cand = sweep(KZ_now[,loc_ind] + Kgam_cand[,rep(1,N_rep)]+ xb,2,
                      alpha_i_now[genus_ind]+alpha_si_now[alp_is_ind_vec] ,"+")
      
      like_dif = sum(dnorm(c(Y),c(mu_cand),mod_var,log = TRUE)) - 
        sum(dnorm(c(Y),c(mu_now),mod_var,log = TRUE)) 
      
      prior_dif = my_dmvnorm(log(v_cand_genus) - m_v_genus_now,v_v_genus_now,
                             log_det_genus,inv_R_genus) - 
        my_dmvnorm(log(v_now_genus) - m_v_genus_now,v_v_genus_now,
                   log_det_genus,inv_R_genus)
      
      if(like_dif + prior_dif + MH_dif > log(runif(1)) &
         !is.na(like_dif + prior_dif + MH_dif)){   
        
        mu_now = mu_cand
        Kgam_now = Kgam_cand
        v_now_genus = v_cand_genus
        K_genus = K_genus_cand
        
        count_v_genus[j] = count_v_genus[j] + 1
      }
    }
  }
  
  v_now = 1/( sum_R_genus/v_v_genus_now + 1/9)
  m_now = sum(R_genus_colsum * log(v_now_genus))/v_v_genus_now + 3/9
  m_v_genus_now = rnorm(1, m_now * v_now, sqrt(v_now))
  
  a_now = 5 + N_knot_genus/2
  b_now = 2 + quad.form(inv_R_genus,log(v_now_genus) - m_v_genus_now)/2
  v_v_genus_now = 1/rgamma(1,a_now,b_now)
  
  v_cand_spat = v_now_spat
  v_cand_spat = rlnorm(1,log(v_now_spat),cand_v_spat)
  
  K_spat_cand = K_spat
  
  MH_dif = dlnorm(v_now_spat,log(v_cand_spat),cand_v_spat,log = TRUE) - 
    dlnorm(v_cand_spat,log(v_now_spat),cand_v_spat,log = TRUE)
  
  K_spat_cand = kern(dif_temp_spat,v_cand_spat)
  KA_cand = (K_spat_cand %*% A_now)
  KZ_cand = KA_cand %*% z_now
  
  mu_cand = sweep(KZ_cand[,loc_ind] + Kgam_now[,rep(1,N_rep)]+ xb,2,
                  alpha_i_now[genus_ind] ,"+")
  
  like_dif = sum(dnorm(c(Y),c(mu_cand),mod_var,log = TRUE)) - 
    sum(dnorm(c(Y),c(mu_now),mod_var,log = TRUE)) 
  
  prior_dif = dgamma(v_cand_spat,5,1/10,log = TRUE) - dgamma(v_now_spat,5,1/10,log = TRUE)
  
  if(like_dif + prior_dif + MH_dif > log(runif(1)) &
     !is.na(like_dif + prior_dif + MH_dif)){
    
    mu_now = mu_cand
    KZ_now = KZ_cand
    KA = KA_cand
    v_now_spat = v_cand_spat
    K_spat = K_spat_cand
    
    count_v_spat = count_v_spat + 1
  }
  
  
  ###################### K_beta
  
  
  v_cand_beta = v_now_beta
  v_cand_beta = rlnorm(1,log(v_now_spat),cand_v_spat)
  
  K_beta_cand = K_beta
  
  MH_dif = dlnorm(v_now_beta,log(v_cand_beta),cand_v_beta,log = TRUE) - 
    dlnorm(v_cand_beta,log(v_now_beta),cand_v_beta,log = TRUE)
  
  K_beta_cand = kern(dif_temp_beta,v_cand_beta)
  xb_cand = K_beta_cand %*% t(X %*%  beta_now)
  mu_cand = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)]+ xb_cand,2,
                  alpha_i_now[genus_ind] ,"+")
  
  like_dif = sum(dnorm(c(Y),c(mu_cand),mod_var,log = TRUE)) - 
    sum(dnorm(c(Y),c(mu_now),mod_var,log = TRUE)) 
  
  prior_dif = dgamma(v_cand_beta,5,1/10,log = TRUE) - dgamma(v_now_spat,5,1/10,log = TRUE)
  
  if(like_dif + prior_dif + MH_dif > log(runif(1)) &
     !is.na(like_dif + prior_dif + MH_dif)){
    
    mu_now = mu_cand
    xb = xb_cand
    v_now_beta = v_cand_beta
    K_beta = K_beta_cand
    
    count_v_beta = count_v_beta + 1
  }
  
  
  
  
  ########   ########   ########   ########   ########   ######## 
  ######## update tau2
  ########   ########   ########   ########   ########   ######## 
  
  for(j in 1:n_par_tau){
    
    beta_tau2_cand = beta_tau2_now
    tau2_cand = tau2_now
    
    beta_tau2_cand[j] = rnorm(1,beta_tau2_now[j],cand_v_tau[j])
    
    tau2_cand = c(exp(K_tau %*% beta_tau2_cand))
    
    if(j == 1){
      prior_dif = dnorm(beta_tau2_cand[j],0,100,log = TRUE) - 
        dnorm(beta_tau2_now[j],0,100,log = TRUE)
    } else{
      prior_dif = dnorm(beta_tau2_cand[j],0,3,log = TRUE) - 
        dnorm(beta_tau2_now[j],0,3,log = TRUE)
    }
    
    like_dif = sum(dnorm(c(Y[ind_tau_spline[[j]],]),c(mu_now[ind_tau_spline[[j]],]),
                         rep(sqrt(tau2_cand[ind_tau_spline[[j]]]),times = N_rep),log = TRUE)) - 
      sum(dnorm(c(Y[ind_tau_spline[[j]],]),c(mu_now[ind_tau_spline[[j]],]),
                rep(sqrt(tau2_now[ind_tau_spline[[j]]]),times = N_rep),log = TRUE)) 
    
    if(like_dif + prior_dif > log(runif(1))){
      
      tau2_now = tau2_cand
      beta_tau2_now = beta_tau2_cand
      
      count_v_tau[j] = count_v_tau[j] + 1
    }
    
    
  }
  
  
  
  
  KtK_spat = quad.form(Matrix::diag(1/tau2_now),K_spat)
  KAtKA_spat = quad.form(Matrix::diag(1/tau2_now),KA)
  KtK_genus= quad.form(Matrix::diag(1/tau2_now),K_genus)
  KtK_beta = quad.form(Matrix::diag(1/tau2_now),K_beta)
  

  
  
  if(i  %% tune == 0){
    
    if(i < burn){
      acc_v_genus = count_v_genus / tune; count_v_genus = rep(0,N_knot_genus)
      acc_v_spat = count_v_spat / tune; count_v_spat = rep(0,1)
      acc_v_beta = count_v_beta / tune; count_v_beta = rep(0,1)
      
      acc_v_tau = count_v_tau / tune; count_v_tau = rep(0,n_par_tau)
      acc_phi = count_phi / tune; count_phi = rep(0,r)
      acc_phi_alp = count_phi_alp / tune; count_phi_alp = 0
      
      cand_v_genus = ifelse(acc_v_genus > 0.6, cand_v_genus*2,
                            ifelse(acc_v_genus < 0.2, cand_v_genus/3, cand_v_genus) )
      cand_v_spat = ifelse(acc_v_spat > 0.6, cand_v_spat*2,
                           ifelse(acc_v_spat < 0.2, cand_v_spat/3, cand_v_spat) )
      cand_v_beta = ifelse(acc_v_beta > 0.6, cand_v_beta*2,
                           ifelse(acc_v_beta < 0.2, cand_v_beta/3, cand_v_beta) )
      cand_v_tau = ifelse(acc_v_tau > 0.6, cand_v_tau*2,
                          ifelse(acc_v_tau < 0.2, cand_v_tau/3, cand_v_tau) )
      
      cand_phi_tau = ifelse(acc_phi > 0.6, cand_phi_tau*2,
                            ifelse(acc_phi < 0.2, cand_phi_tau/3, cand_phi_tau) )
      
      cand_phi_tau_alp = ifelse(acc_phi_alp > 0.6, cand_phi_tau_alp*2,
                            ifelse(acc_phi_alp < 0.2, cand_phi_tau_alp/3, 
                                   cand_phi_tau_alp) )
      
    } 
    
    time_its <- (proc.time() - st)[3] / (i - 1)
    time_used <- round((proc.time() - st)[3]/(60),digits=4)
    time_left <- round(time_its * (reps*thin + burn - i )/(60),digits=4)
    cat("\r", i, " of ", reps*thin + burn,"||| Time left: ",floor(time_left/60),
        " hours",time_left%%60," minutes")# |||| like = ", log_lik[i-1]) 
    flush.console()
    
    
  }
  

  
  
  if(i > burn & i %% thin == 0){
    ### save parameters
    
    
    alpha_si[(i - burn)/thin,] = alpha_si_now
    sig2_alp_s[(i - burn)/thin] = sig2_alp_s_now
    
    alpha_i[(i - burn)/thin,] = alpha_i_now
    alpha[(i - burn)/thin] = alpha_now
    beta[(i - burn)/thin,,] = beta_now
    gamma[(i - burn)/thin,] = gamma_now
    Sigma_gam[(i - burn)/thin,,] = Sigma_gam_now
    A[(i - burn)/thin,,] = A_now
    z[(i - burn)/thin,,] = z_now
    phi_spat[(i - burn)/thin,] = phi_spat_now
    phi_alp_s[(i - burn)/thin] = phi_alp_s_now
    
    beta_tau2[(i - burn)/thin,] = beta_tau2_now
    sig2_alp[(i - burn)/thin] = sig2_alp_now
    sig2_gamma[(i - burn)/thin] = sig2_gamma_now
    sig2_A[(i-burn)/thin] = sig2_A_now
    
    v_genus[(i - burn)/thin,] = v_now_genus 
    v_spat[(i - burn)/thin] = v_now_spat 
    v_beta[(i - burn)/thin] = v_now_beta 
    
    m_v_genus[(i - burn)/thin] = m_v_genus_now
    v_v_genus[(i - burn)/thin] = v_v_genus_now

    
    mod_var = rep(sqrt(tau2_now),times = N_rep)
    like_save[(i - burn)/thin] =  sum(dnorm(c(Y),c(mu_now),mod_var,log = TRUE))
  }
  
  

}


save.image("aster.RData")

