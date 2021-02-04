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
library(viridis)
rm(list = ls())
############## map for families


try(setwd("C:/Users/philaw/Box/Research/ERT/code/test_code/final2"),silent = TRUE)

# make_the_coef_plot = function(beta_funcs_mean,CI, qq,col_use){
#   plot(450:949,beta_funcs_mean[,qq],ylab = paste("Effect of",colnames(X)[qq]),
#        xlab = "wavelength (nm)",cex.lab = 1.6,ylim = range(CI[,,qq]),col = "white",
#        xaxs="i", yaxs="i",cex.axis = 1.6)
#   polygon(c(450:949,949:450), c(CI[1,,qq],CI[2,N_wave:1,qq]),
#           col = "lightgray",border = col_use,lwd = 2)
#   lines(450:949,beta_funcs_mean[,qq],col = col_use,lwd = 2)
# }
# 
# 
# 
# make_the_coef_lines = function(beta_funcs_mean,CI, qq){
#   lines(450:949,beta_funcs_mean[,qq],ylab = paste("Effect of",colnames(X)[qq]),
#        xlab = "wavelength (nm)",cex.lab = 1.6,ylim = range(CI[,,qq]),col = "white",
#        xaxs="i", yaxs="i",cex.axis = 1.6)
#   polygon(c(450:949,949:450), c(CI[1,,qq],CI[2,N_wave:1,qq]),
#           col = "lightgray",border = "red",lwd = 2)
#   lines(450:949,beta_funcs_mean[,qq],col = "cornflowerblue",lwd = 2)
# }


quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

N_wave = 500


load("aster.RData")

attach(stan_data)


f_X = X[loc_ind,] %*% solve( t(X[loc_ind,]) %*% X[loc_ind,])
f_K = solve( t(K_beta) %*% K_beta) %*% t(K_beta)
########

bet_funcs_aster_confound = lapply(1:reps,function(x){K_beta %*% t(beta[x,,])})


# proj_X = f_X %*% t(X[loc_ind,])
# proj_K =  K_beta %*% f_K
# 
# 
# 
# idx_use = which(1:reps  %% 10 == 0)
# 
# prop_confound = matrix(0,length(idx_use),3)  ##### alpha_i(s), gamma(t), eta(s,t)
# var_left_orth = matrix(0,length(idx_use),5)  ##### unexplained, beta, alpha_i(s), gamma(t), eta(s,t)
# var_left_raw = matrix(0,length(idx_use),5)  ##### unexplained, beta, alpha_i(s), gamma(t), eta(s,t)
# 
# colnames(var_left_orth) <- colnames(var_left_raw) <- c("all terms","no beta","no genus","no gamma","no eta")
# 
# var_wave = array(0,c(length(idx_use),5,N_wave))
# var_site = array(0,c(length(idx_use),5,N_rep))
# 
# 
# 
# 
# for(i in 1:length(idx_use)){
#   
#   for(j in 1:N_knot_genus){
#     K_genus[,j] = kern(dif_temp_genus[,j],v_genus[idx_use[i],j])
#   }
#   
#   K_spat = kern(dif_temp_spat,v_spat[idx_use[i]])
#   KA = (K_spat %*% A[idx_use[i],,])
#   
#   ### Raw quantities
#   
#   Kgam_now = K_genus %*% gamma[idx_use[i],]
#   
#   KZ_now = KA %*% z[idx_use[i],,]
#   
#   xb_now = K_beta %*% t(beta[idx_use[i],,]) %*% t(X[loc_ind,])
#   
#   ### Projections / confounding
#   
#   alp_raw = sweep(matrix(0,N_wave,N_rep),2,
#                   alpha_i[idx_use[i],genus_ind] + alpha_si[idx_use[i],alp_is_ind_vec] - alpha[idx_use[i]],
#                   "+")
#   
#   rand_eff = KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)] + alp_raw
#   
#   beta_update = K_beta %*% (f_K %*% (rand_eff %*% f_X))
#   rand_eff_proj = beta_update %*% t(X[loc_ind,])
#   
#   xb_proj = xb_now + rand_eff_proj
#   
#   KZ_proj = proj_K %*% (KZ_now[,loc_ind] %*% proj_X)
#   Kgam_proj = proj_K %*% (Kgam_now[,rep(1,N_rep)] %*% proj_X)
#   genus_stuff =  proj_K %*% (alp_raw %*% proj_X)
#   
#   # prop_confound[i,] = c(var(c(genus_stuff))/var(c(alp_raw)),
#   #                       var(c(Kgam_proj))/var(c(Kgam_now)),
#   #                       var(c(KZ_proj))/var(c(KZ_now)))
#   
#   genus_orth =  alp_raw - genus_stuff
#   KZ_orth = KZ_now[,loc_ind] - KZ_proj
#   Kgam_orth = Kgam_now[,rep(1,N_rep)] - Kgam_proj
#   
#   #### Proportion of Proportion Variance Explained
#   
#   mu_now_all1 = alpha[idx_use[i]] + xb_proj + genus_orth + Kgam_orth + KZ_orth
#   # mu_now_beta = alpha[idx_use[i]]  + genus_orth + Kgam_orth + KZ_orth
#   # mu_now_genus = alpha[idx_use[i]] + xb_proj  + Kgam_orth + KZ_orth
#   # mu_now_gamma = alpha[idx_use[i]] + xb_proj + genus_orth  + KZ_orth
#   # mu_now_eta = alpha[idx_use[i]] + xb_proj + genus_orth + Kgam_orth
#   
#   var_left_orth[i,] = c(var(c(Y - mu_now_all1)),
#                         var(c(xb_proj)),
#                         var(c(genus_orth)),
#                         var(c(Kgam_orth)),
#                         var(c(KZ_orth)))
#   
#   var_wave[i,,] = rbind(apply(Y - mu_now_all1,1,var),
#                         apply(xb_proj,1,var),
#                         apply(genus_orth,1,var),
#                         apply(Kgam_orth,1,var),
#                         apply(KZ_orth,1,var))
#   
#   
#   var_site[i,,] = rbind(apply(Y - mu_now_all1,2,var),
#                         apply(xb_proj,2,var),
#                         apply(genus_orth,2,var),
#                         apply(Kgam_orth,2,var),
#                         apply(KZ_orth,2,var))
#   
#   mu_now_all2 = alpha[idx_use[i]] + xb_now + alp_raw + Kgam_now[,rep(1,N_rep)] + KZ_now[,loc_ind]
#   # mu_now_beta = alpha[idx_use[i]] + alp_raw + Kgam_now[,rep(1,N_rep)] + KZ_now[,loc_ind]
#   # mu_now_genus = alpha[idx_use[i]] + xb_now + Kgam_now[,rep(1,N_rep)] + KZ_now[,loc_ind]
#   # mu_now_gamma = alpha[idx_use[i]] + xb_now + alp_raw  + KZ_now[,loc_ind]
#   # mu_now_eta = alpha[idx_use[i]] + xb_now + alp_raw + Kgam_now[,rep(1,N_rep)]
#   
#   var_left_raw[i,] = c(var(c(Y - mu_now_all2)),
#                        var(c(xb_now)),
#                        var(c(alp_raw)),
#                        var(c(Kgam_now)),
#                        var(c( KZ_now[,loc_ind])))
#   
#   
#   
#   cat("\r", i, " of ", length(idx_use))
#   
# }
# 
# var_left_orth_scale = var_left_orth
# var_left_raw_scale = var_left_raw
# 
# for(i in 1:length(idx_use)){
#   
#   var_left_raw_scale[i,] = var_left_raw[i,] / sum(var_left_raw[i,])
#   var_left_orth_scale[i,1] = var_left_raw_scale[i,1]
#   var_left_orth_scale[i,-1] = (1 - var_left_orth_scale[i,1]) * var_left_orth[i,-1] / sum(var_left_orth[i,-1]) 
#   
# }
# 
# var_summaries = list(var_left_orth = var_left_orth_scale,
#                      var_left_raw = var_left_raw_scale,
#                      var_wave= simplify2array(lapply(1:length(idx_use),function(x){
#                        sweep(var_wave[x,,],2,apply(var_wave[x,,],2,sum),"/")
#                      })))
# saveRDS(var_summaries,"var_summaries.rds")

# colnames(prop_confound) = c("genus","wavelength","space-wavelength")
# df_temp = melt(prop_confound)
# colnames(df_temp) = c("rep","Random Effect","Proportion Confounded")
# 
# pdf("prop_confound.pdf")
# ggplot(data = df_temp) +
# geom_boxplot(aes(x = `Random Effect`,y = `Proportion Confounded`,col = `Random Effect`,fill = `Random Effect`)) +
#   theme(axis.text.x = element_text(angle = 325, vjust = 1, hjust=0))+
#   theme(axis.text=element_text(size=24),
#         # axis.text.x = element_blank(),
#         axis.title=element_text(size=24),
#         legend.position = "none") +# ylim(c(0,.25)) +
#   scale_x_discrete(labels = c('genus' = expression(paste(alpha[i],"(s) - ",alpha)),
#                               'wavelength' = expression(paste(gamma,"(t)")),
#                               'space-wavelength' = expression(paste(eta,"(s,t)"))))
# 
# dev.off()

var_summaries = readRDS("var_summaries.rds")
attach(var_summaries)

apply(var_left_orth,2,mean)
apply(var_left_orth,2,sd)
apply(var_left_orth,2,quantile,c(.025,.975))

var_left_orth_trim  = apply(var_left_orth,2,function(x){
  qq = quantile(x,c(.025,.975))
  x[which(x <= qq[2]   & x >= qq[1])]
})

var_left_raw_trim  = apply(var_left_raw,2,function(x){
  qq = quantile(x,c(.025,.975))
  x[which(x <= qq[2]   & x >= qq[1])]
})

df_temp1 = as.data.frame(var_left_orth_trim)
colnames(df_temp1) = c("error","Covariates","genus","gamma(t)","eta(s,t)")
df_temp1 = reshape::melt(df_temp1)

df_temp2 = as.data.frame(var_left_raw_trim)
colnames(df_temp2) = c("error","Covariates","genus","gamma(t)","eta(s,t)")
df_temp2 = reshape::melt(df_temp2)

df_temp = rbind(df_temp1,df_temp2)
df_temp$orth = rep(c("orthogonalized","confounded"),each = 5* nrow(var_left_raw_trim))

colnames(df_temp) = c("model term","Proportion Variance Explained","Orthogonalized")


pdf("var_explained.pdf",width =10,height = 6)
ggplot(data = df_temp) + 
  geom_boxplot(aes(x = `model term`,y = `Proportion Variance Explained`,
                   col = `Orthogonalized`,fill = `Orthogonalized`)) +
  theme(axis.text.x = element_text(angle = 325, vjust = 1, hjust=0))+
  theme(axis.text=element_text(size=24),
        # axis.text.x = element_blank(),
        axis.title=element_text(size=24),
        legend.title=element_text(size=22),
        legend.text=element_text(size=22)) + 
  scale_x_discrete(labels = c('error' = expression(paste(epsilon[ij],"(s,t)")),
                              'Covariates' = expression(paste("x(s)",beta,"(t)")),
                              'genus' = expression(paste(alpha[i],"(s) + ",alpha[i]," - ",alpha)),
                              'gamma(t)' = expression(paste(gamma,"(t)")),
                              'eta(s,t)' = expression(paste(eta,"(s,t)"))))
dev.off()

col_use = viridis(5)
library(reshape2)

dimnames(var_wave) = list(c("error","covariates","genus","gamma(t)","eta(s,t)"),
                          450:949,1:dim(var_wave)[3])

# for(i in 1:5){
#   for(j in 1:500){
#     var_wave[i,j,] = var_wave[i,j,]/(1 - var_wave[4,j,])
#   }
# }
# 
# var_wave = var_wave[-4,,]

temp1 = melt(apply(var_wave,c(1,2),mean))
temp2 = melt(apply(var_wave,c(1,2),quantile,0.025))
temp3 = melt(apply(var_wave,c(1,2),quantile,0.975))
temp = rbind(temp1,temp2,temp3)
temp$line = rep(c("mean","2.5%","97.5%"),each = 2500)

colnames(temp) = c("model term","Wavelength","Proportion Variance Explained","line")
temp$group = factor(temp$`line`):temp$`model term`
# plot(temp[1,],type = "l",lwd =3,lty = 1,ylim = range(temp),col = )     ### error
# lines(temp[2,],type = "l",lwd =3,lty = 2)                  ### xb  
# lines(temp[3,],type = "l",lwd =3,lty = 3)                  ### genus
# lines(temp[4,],type = "l",lwd =3,lty = 4)                  ### gamma(t)
# lines(temp[5,],type = "l",lwd =3,lty =5)                   ### eta(s,t)
# 

temp_CI = temp[temp$line != "mean",]
temp_CI = temp_CI[c(order(temp$Wavelength[temp$line == "2.5%"]),
                    2500 + order(temp$Wavelength[temp$line == "97.5%"],decreasing = TRUE)),]

pdf("between_site_wavelength.pdf",width =12, height=6)
ggplot() + 
  geom_polygon(data = temp_CI,
               aes(x = Wavelength,y = `Proportion Variance Explained`,
                   col = `model term`,group = `model term`),
               size=1.3,alpha=.2) + 
  geom_line(data = temp[temp$line == "mean",],
            aes(x = Wavelength,y = `Proportion Variance Explained`,
                col = `model term`,group = `model term`),
            size=1.3,alpha=.7) +
  #theme(axis.text.x = element_text(angle = 325, vjust = 1, hjust=0))+
  theme(axis.text=element_text(size=24),
        # axis.text.x = element_blank(),
        axis.title=element_text(size=24),
        legend.title=element_text(size=22),
        legend.text=element_text(size=22)) + 
  scale_color_viridis_d(labels = c('error' = expression(paste(epsilon[ij],"(s,t)")),
                                   'covariates' = expression(paste("x(s)",beta,"(t)")),
                                   'genus' = expression(paste(alpha[i],"(s) + ",alpha[i]," - ",alpha)),
                                   'gamma(t)' = expression(paste(gamma,"(t)")),
                                   'eta(s,t)' = expression(paste(eta,"(s,t)"))),
                        alpha = 0.7) + 
  theme(legend.text.align = 0) +
  scale_x_continuous(expand = c(0, 0)) 
dev.off()

#beta_spat_un <- 
beta_both_un <- 
  bet_funcs_aster_confound

tau2_aster_save = matrix(0,reps,N_wave)


for(i in 1:reps){
  
  for(j in 1:N_knot_genus){
    K_genus[,j] = kern(dif_temp_genus[,j],v_genus[i,j])
  }
  
  K_spat = kern(dif_temp_spat,v_spat[i])
  
  Kgam_now = K_genus %*% gamma[i,]
  
  KA = (K_spat %*% A[i,,])
  KZ_now = KA %*% z[i,,]
  
  
  
  eta_now = sweep(KZ_now[,loc_ind] + Kgam_now[,rep(1,N_rep)],2,
                  alpha_i[i,genus_ind] - alpha[i] + alpha_si[i,alp_is_ind_vec] ,"+")
  
  
  beta_both_un[[i]] = bet_funcs_aster_confound[[i]] + K_beta %*% (f_K %*% (eta_now %*% f_X))
  # beta_spat_un[[i]] = bet_funcs_aster_confound[[i]] + (eta_now %*% f_X)
  
  tau2_aster_save[i,]  = c(exp(K_tau %*% beta_tau2[i,]))
  
  
  cat("\r", i, " of ", reps) 
  
}

lm_temp = lm(c(Y) ~ X[loc_ind,] %x% K_beta)

bet_coef_ls = K_beta %*% matrix(coef(lm_temp)[-1],nrow =ncol(K_beta))


plot(bet_coef_ls[,1])

################# look at beta functions

colnames(X) = c("Elevation","Annual Precipitation","Avg Min Jan Temp", "Rainfall Concentration")

bet_funcs_mean_aster_confound = Reduce("+", bet_funcs_aster_confound) / length(bet_funcs_aster_confound)
CI_aster_confound = apply(simplify2array(bet_funcs_aster_confound), 1:2, quantile, prob = c(0.025, 0.975))


bet_funcs_mean_aster_both = Reduce("+", beta_both_un) / length(beta_both_un)
CI_aster_un = apply(simplify2array(beta_both_un), 1:2, quantile, prob = c(0.025, 0.975))

# bet_funcs_mean_aster_spat = Reduce("+", beta_spat_un) / length(beta_spat_un)
# CI_aster_spat = apply(simplify2array(beta_spat_un), 1:2, quantile, prob = c(0.025, 0.975))


pdf("elevation.pdf",width = 12,height = 8)

par(mar = c(5,5,1,1))
qq = 1
col_use = "black"

plot(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
     xlab = "wavelength (nm)",cex.lab = 2.4,ylim = range(CI_aster_un[,,qq]),col = "white",
     xaxs="i", yaxs="i",cex.axis = 1.6)
lines(450:949,bet_coef_ls[,qq],col = "white",lwd = 2)



col_use = "red"

lines(450:949,bet_funcs_mean_aster_both[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = "white",
      xaxs="i", yaxs="i",cex.axis = 1.6)
polygon(c(450:949,949:450), c(CI_aster_un[1,,qq],CI_aster_un[2,N_wave:1,qq]),
        col = "white",border = col_use,lwd = 3,lty = 2)
lines(450:949,bet_funcs_mean_aster_both[,qq],col = "gray",lwd = 3)


col_use = "black"

lines(450:949,bet_coef_ls[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = col_use,lwd = 3,lty = 4)

abline(h = 0,col = "blue",lty = 3, lwd = 3)

legend("topleft",c("Posterior Mean","95% Credible Interval",expression(paste(beta,"(t) without Random Effects"))),
       col = c("gray","red","black"),lty = c(1,2,4),lwd = 4,cex = 2)


dev.off()

pdf("precip.pdf",width = 12,height = 8)
par(mar = c(5,5,1,1))

qq = 2
col_use = "black"

plot(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
     xlab = "wavelength (nm)",cex.lab = 2.4,ylim = range(CI_aster_un[,,qq]),col = "white",
     xaxs="i", yaxs="i",cex.axis = 1.6)
lines(450:949,bet_coef_ls[,qq],col = "white",lwd = 2)



col_use = "red"

lines(450:949,bet_funcs_mean_aster_both[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = "white",
      xaxs="i", yaxs="i",cex.axis = 1.6)
polygon(c(450:949,949:450), c(CI_aster_un[1,,qq],CI_aster_un[2,N_wave:1,qq]),
        col = "white",border = col_use,lwd = 3,lty = 2)
lines(450:949,bet_funcs_mean_aster_both[,qq],col = "gray",lwd = 3)


col_use = "black"

lines(450:949,bet_coef_ls[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = col_use,lwd = 4,lty = 4)

abline(h = 0,col = "blue",lty = 3, lwd = 3)

legend("topleft",c("Posterior Mean","95% Credible Interval",expression(paste(beta,"(t) without Random Effects"))),
       col = c("gray","red","black"),lty = c(1,2,4),lwd = 4,cex = 2)

dev.off()

pdf("temp.pdf",width = 12,height = 8)
par(mar = c(5,5,1,1))

qq = 3

col_use = "black"

plot(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
     xlab = "wavelength (nm)",cex.lab = 2.4,ylim = range(CI_aster_un[,,qq]),col = "white",
     xaxs="i", yaxs="i",cex.axis = 1.6)
lines(450:949,bet_coef_ls[,qq],col = "white",lwd = 2)



col_use = "red"

lines(450:949,bet_funcs_mean_aster_both[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = "white",
      xaxs="i", yaxs="i",cex.axis = 1.6)
polygon(c(450:949,949:450), c(CI_aster_un[1,,qq],CI_aster_un[2,N_wave:1,qq]),
        col = "white",border = col_use,lwd = 3,lty = 2)
lines(450:949,bet_funcs_mean_aster_both[,qq],col = "gray",lwd = 3)


col_use = "black"

lines(450:949,bet_coef_ls[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = col_use,lwd = 4,lty = 4)

abline(h = 0,col = "blue",lty = 3, lwd = 3)

legend("topleft",c("Posterior Mean","95% Credible Interval",expression(paste(beta,"(t) without Random Effects"))),
       col = c("gray","red","black"),lty = c(1,2,4),lwd = 4,cex = 2)
dev.off()

pdf("rainfallconcentration.pdf",width = 12,height = 8)
par(mar = c(5,5,1,1))

qq = 4

col_use = "black"

plot(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
     xlab = "wavelength (nm)",cex.lab = 2.4,ylim = range(CI_aster_un[,,qq]),col = "white",
     xaxs="i", yaxs="i",cex.axis = 1.6)
lines(450:949,bet_coef_ls[,qq],col = "white",lwd = 2)



col_use = "red"

lines(450:949,bet_funcs_mean_aster_both[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = "white",
      xaxs="i", yaxs="i",cex.axis = 1.6)
polygon(c(450:949,949:450), c(CI_aster_un[1,,qq],CI_aster_un[2,N_wave:1,qq]),
        col = "white",border = col_use,lwd = 3,lty = 2)
lines(450:949,bet_funcs_mean_aster_both[,qq],col = "gray",lwd = 3)


col_use = "black"

lines(450:949,bet_coef_ls[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = col_use,lwd = 4,lty = 4)

abline(h = 0,col = "blue",lty = 3, lwd = 3)
legend("bottomleft",c("Posterior Mean","95% Credible Interval",expression(paste(beta,"(t) without Random Effects"))),
       col = c("gray","red","black"),lty = c(1,2,4),lwd = 4,cex = 2)
dev.off()






pdf("elevation_confound.pdf",width = 12,height = 8)
par(mar = c(5,5,1,1))

qq = 1

col_use = "black"

plot(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
     xlab = "wavelength (nm)",cex.lab = 2.4,ylim = range(CI_aster_confound[,,qq]),col = "white",
     xaxs="i", yaxs="i",cex.axis = 1.6)
lines(450:949,bet_coef_ls[,qq],col = "white",lwd = 2)



col_use = "red"

lines(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = "white",
      xaxs="i", yaxs="i",cex.axis = 1.6)
polygon(c(450:949,949:450), c(CI_aster_confound[1,,qq],CI_aster_confound[2,N_wave:1,qq]),
        col = "white",border = col_use,lwd = 3,lty = 2)
lines(450:949,bet_funcs_mean_aster_confound[,qq],col = "gray",lwd = 3)

abline(h =0,lwd = 3,lty = 3)

dev.off()

pdf("precip_confound.pdf",width = 12,height = 8)
par(mar = c(5,5,1,1))

qq = 2

col_use = "black"

plot(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
     xlab = "wavelength (nm)",cex.lab = 2.4,ylim = range(CI_aster_confound[,,qq]),col = "white",
     xaxs="i", yaxs="i",cex.axis = 1.6)
lines(450:949,bet_coef_ls[,qq],col = "white",lwd = 2)



col_use = "red"

lines(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = "white",
      xaxs="i", yaxs="i",cex.axis = 1.6)
polygon(c(450:949,949:450), c(CI_aster_confound[1,,qq],CI_aster_confound[2,N_wave:1,qq]),
        col = "white",border = col_use,lwd = 3,lty = 2)
lines(450:949,bet_funcs_mean_aster_confound[,qq],col = "gray",lwd = 3)

abline(h =0,lwd = 3,lty = 3)

dev.off()

pdf("temp_confound.pdf",width = 12,height = 8)
par(mar = c(5,5,1,1))

qq = 3

col_use = "black"

plot(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
     xlab = "wavelength (nm)",cex.lab = 2.4,ylim = range(CI_aster_confound[,,qq]),col = "white",
     xaxs="i", yaxs="i",cex.axis = 1.6)
lines(450:949,bet_coef_ls[,qq],col = "white",lwd = 2)



col_use = "red"

lines(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = "white",
      xaxs="i", yaxs="i",cex.axis = 1.6)
polygon(c(450:949,949:450), c(CI_aster_confound[1,,qq],CI_aster_confound[2,N_wave:1,qq]),
        col = "white",border = col_use,lwd = 3,lty = 2)
lines(450:949,bet_funcs_mean_aster_confound[,qq],col = "gray",lwd = 3)

abline(h =0,lwd = 3,lty = 3)

dev.off()


pdf("rainfallconcentration_confound.pdf",width = 12,height = 8)
par(mar = c(5,5,1,1))

qq = 4

col_use = "black"

plot(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
     xlab = "wavelength (nm)",cex.lab = 2.4,ylim = range(CI_aster_confound[,,qq]),col = "white",
     xaxs="i", yaxs="i",cex.axis = 1.6)
lines(450:949,bet_coef_ls[,qq],col = "white",lwd = 2)



col_use = "red"

lines(450:949,bet_funcs_mean_aster_confound[,qq],ylab = paste("Effect of",colnames(X)[qq]),
      xlab = "wavelength (nm)",cex.lab = 1.6,col = "white",
      xaxs="i", yaxs="i",cex.axis = 1.6)
polygon(c(450:949,949:450), c(CI_aster_confound[1,,qq],CI_aster_confound[2,N_wave:1,qq]),
        col = "white",border = col_use,lwd = 3,lty = 2)
lines(450:949,bet_funcs_mean_aster_confound[,qq],col = "gray",lwd = 3)

abline(h =0,lwd = 3,lty = 3)

dev.off()


##################################################################################
#########################################
##################################################################################

temp = t(apply(simplify2array(beta_both_un),c(2,3),function(x){mean(abs(x))}))
colnames(temp) = colnames(X)
imp_avg = apply(temp,2,mean)


xtable(t(apply(temp,2,function(x){c(mean(x),quantile(x,c(.025,.975)))})),digits =3)

var_imp = melt(temp)
colnames(var_imp) = c("reps","covariate","Relative Importance")


# ggplot(var_imp, aes(x=covariate, y=`Relative Importance`, fill=covariate)) +
#   guides(fill=F) +
#   coord_flip() +
#   stat_summary(fun.data = quantiles_95, geom="boxplot")+
#   # theme(axis.text.x = element_text(angle = 270, vjust = 0.2, hjust=0))+
#   theme(axis.text=element_text(size=16),
#         axis.title=element_text(size=18),
#         legend.title=element_text(size=16),
#         legend.text=element_text(size=14))

pdf("aster_var_imp.pdf")
ggplot(data = var_imp) + 
  geom_boxplot(aes(x = covariate,y = `Relative Importance`,col = covariate,fill = covariate)) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.2, hjust=0))+
  theme(axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        axis.title=element_text(size=18),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16))
dev.off()


temp = t(apply(simplify2array(beta_both_un),c(2,3),function(x){mean(x)}))
colnames(temp) = colnames(X)
bet_avg = apply(temp,2,mean)
xtable(t(apply(temp,2,function(x){c(mean(x),quantile(x,c(.025,.975)))})),digits =3)


var_avg = melt(temp)
colnames(var_avg) = c("reps","covariate","Average Coefficient")

pdf("aster_var_avg.pdf")
ggplot(data = var_avg) + 
  geom_boxplot(aes(x = covariate,y = `Average Coefficient`,
                   col = covariate,fill = covariate)) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.2, hjust=0))+
  theme(axis.text.y=element_text(size=16),
        axis.text.x = element_blank(),
        axis.title=element_text(size=18),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16))
dev.off()

# ggplot(var_avg, aes(x=covariate, y=`Average Coefficient`, fill=covariate)) +
#   guides(fill=F) +
#   coord_flip() +
#   stat_summary(fun.data = quantiles_95, geom="boxplot")+
#   # theme(axis.text.x = element_text(angle = 270, vjust = 0.2, hjust=0))+
#   theme(axis.text=element_text(size=16),
#         axis.title=element_text(size=18),
#         legend.title=element_text(size=16),
#         legend.text=element_text(size=14))





prop_sig = matrix(apply(CI_aster_un,c(3),function(x){
  mean(apply(x,2,function(y){
    1 *( y[2] < 0 | y[1] > 0 )
  }))
}),nrow = 1)

colnames(prop_sig) = colnames(X)
prop_sig 
bet_sum = cbind(t(prop_sig),imp_avg,bet_avg)


xtable(bet_sum,digits = 4)
