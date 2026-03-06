library(ggplot2)
library(dplyr)

# scenario 'correctModel' or 'negbinomModel'
scenario = 'correctModel'

# load T = 30
file_30 = paste("~/Downloads/SimStudies_with_RData/SimStudy_",scenario,"/",
           "SimStudy_",scenario,"_T30/SimStudy_",scenario,"_T30.RData",sep='')
load(file_30)

times = c()
case_time = c()

case = c()
point_est = c()
param = c()
model = c()
CI95_inf = c()
CI95_sup = c()

eff_size = c()

# Case 1
for (i in 1:length(res_case1)){
  times = c(times,res_case1[[i]]$compTimes)
  case_time = c(case_time, rep(1,length(res_case1[[i]]$compTimes)))
  
  CI95_inf = c(CI95_inf,res_case1[[i]]$CI95[,1])
  CI95_sup = c(CI95_sup,res_case1[[i]]$CI95[,2])
  point_est = c(point_est ,res_case1[[i]]$pointEst)
  param = c(param,rep(c('theta1','theta2','b'), each =3))
  model = c(model, rep(c('mle','Gibbs','Stan'), 3))
  # ESS
  vec = c(res_case1[[i]]$effSize)
  eff_size = c(eff_size, c(NA,vec[1:2],NA,vec[3:4],NA,vec[5:6]))
  case = c(case,rep(1,length(res_case1[[i]]$pointEst)))
}

# Case 2
for (i in 1:length(res_case2)){
  times = c(times,res_case2[[i]]$compTimes)
  case_time = c(case_time, rep(2,length(res_case2[[i]]$compTimes)))
  
  CI95_inf = c(CI95_inf,res_case2[[i]]$CI95[,1])
  CI95_sup = c(CI95_sup,res_case2[[i]]$CI95[,2])
  point_est = c(point_est ,res_case2[[i]]$pointEst)
  param = c(param,rep(c('theta1','theta2','b'), each =3))
  model = c(model, rep(c('mle','Gibbs','Stan'), 3))
  # ESS
  vec = c(res_case2[[i]]$effSize)
  eff_size = c(eff_size, c(NA,vec[1:2],NA,vec[3:4],NA,vec[5:6]))
  case = c(case,rep(2,length(res_case2[[i]]$pointEst)))
}

times_df = data.frame(CompTimes = times, model = names(times),
                      case = case_time)
points_df = data.frame(PointEst = point_est, param = param, model = model,
                       EffSize = eff_size, case = case,
                       CI95_inf = CI95_inf ,CI95_sup = CI95_sup)
points_df$real = numeric(length(points_df$param))

points_df$real[points_df$param == 'theta1' & case ==1] = theta_case1[1]
points_df$real[points_df$param == 'theta2'& case ==1] = theta_case1[2]
points_df$real[points_df$param == 'b'& case ==1] = theta_case1[3]

points_df$real[points_df$param == 'theta1' & case ==2] = theta_case2[1]
points_df$real[points_df$param == 'theta2'& case ==2] = theta_case2[2]
points_df$real[points_df$param == 'b'& case ==2] = theta_case2[3]

times_df$nsample = rep(30,dim(times_df)[1])
points_df$nsample = rep(30,dim(points_df)[1])

# load T = 100 
file = paste("~/Downloads/SimStudies_with_RData/SimStudy_",scenario,"/",
             "SimStudy_",scenario,"_T100/SimStudy_",scenario,"_T100.RData",
             sep='')
load(file)

times = c()
case_time = c()

case = c()
point_est = c()
param = c()
model = c()
CI95_inf = c()
CI95_sup = c()
eff_size = c()

# Case 1
for (i in 1:length(res_case1)){
  times = c(times,res_case1[[i]]$compTimes)
  case_time = c(case_time, rep(1,length(res_case1[[i]]$compTimes)))
  
  CI95_inf = c(CI95_inf,res_case1[[i]]$CI95[,1])
  CI95_sup = c(CI95_sup,res_case1[[i]]$CI95[,2])
  point_est = c(point_est ,res_case1[[i]]$pointEst)
  param = c(param,rep(c('theta1','theta2','b'), each =3))
  model = c(model, rep(c('mle','Gibbs','Stan'), 3))
  # ESS
  vec = c(res_case1[[i]]$effSize)
  eff_size = c(eff_size, c(NA,vec[1:2],NA,vec[3:4],NA,vec[5:6]))
  case = c(case,rep(1,length(res_case1[[i]]$pointEst)))
}

# Case 2
for (i in 1:length(res_case2)){
  times = c(times,res_case2[[i]]$compTimes)
  case_time = c(case_time, rep(2,length(res_case2[[i]]$compTimes)))
  
  CI95_inf = c(CI95_inf,res_case2[[i]]$CI95[,1])
  CI95_sup = c(CI95_sup,res_case2[[i]]$CI95[,2])
  point_est = c(point_est ,res_case2[[i]]$pointEst)
  param = c(param,rep(c('theta1','theta2','b'), each =3))
  model = c(model, rep(c('mle','Gibbs','Stan'), 3))
  # ESS
  vec = c(res_case2[[i]]$effSize)
  eff_size = c(eff_size, c(NA,vec[1:2],NA,vec[3:4],NA,vec[5:6]))
  case = c(case,rep(2,length(res_case2[[i]]$pointEst)))
}

times_df_100 = data.frame(CompTimes = times, model = names(times),
                      case = case_time)
points_df_100 = data.frame(PointEst = point_est, param = param, model = model,
                       EffSize = eff_size, case = case,
                       CI95_inf = CI95_inf ,CI95_sup = CI95_sup)
points_df_100$real = numeric(length(points_df_100$param))

points_df_100$real[points_df_100$param == 'theta1' & case ==1] = theta_case1[1]
points_df_100$real[points_df_100$param == 'theta2'& case ==1] = theta_case1[2]
points_df_100$real[points_df_100$param == 'b'& case ==1] = theta_case1[3]

points_df_100$real[points_df_100$param == 'theta1' & case ==2] = theta_case2[1]
points_df_100$real[points_df_100$param == 'theta2'& case ==2] = theta_case2[2]
points_df_100$real[points_df_100$param == 'b'& case ==2] = theta_case2[3]


times_df_100$nsample = rep(100,dim(times_df_100)[1])
points_df_100$nsample = rep(100,dim(points_df_100)[1])

all_times_df = rbind(times_df,times_df_100)
all_times_df$case[all_times_df$case==1] = 'High'
all_times_df$case[all_times_df$case==2] = 'Moderate'
all_times_df$case = as.factor(all_times_df$case)
all_times_dfnsample = as.factor(all_times_df$nsample)

all_points_df = rbind(points_df,points_df_100)
all_points_df$case[all_points_df$case==1] = 'High'
all_points_df$case[all_points_df$case==2] = 'Moderate'
all_points_df$case = as.factor(all_points_df$case)
all_points_df$nsample = as.factor(all_points_df$nsample)

rm(times_df,times_df_100,points_df,points_df_100)
###################################### Plots ##################################

#  CI95
all_points_df$isin = as.integer(all_points_df$real>= all_points_df$CI95_inf &
                                  all_points_df$real<=all_points_df$CI95_sup)

summ = all_points_df %>% group_by(param,model,case,nsample)  %>% 
  summarise(cov = sum(isin)/500)
#summ = summ %>% filter(model!='Stan')
levels(summ$case) = c('High', 'Moderate')

summ$nsample = as.factor(summ$nsample)
summ$Model = as.factor(summ$model)

coverage_plot = ggplot(summ) +
  geom_point(aes(x = case, y=cov, color = Model, shape = Model),alpha =0.9) + 
  theme_bw() + facet_grid(nsample ~param ) + geom_hline(aes(yintercept = 0.95),
                                       linetype = "dashed", color = "grey") + 
  xlab('Correlation scenario') + ylab('Coverage')  +   
  scale_color_manual(values = c("brown1","deepskyblue","chartreuse3"))+ 
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position = 'bottom')

ggsave(coverage_plot,file = paste0('./Coverage_plot_',scenario,'.png'),
       width = 6, height = 7)

# MSE of the point estimates
all_points_df = all_points_df %>% mutate(sq_error = (PointEst - real)^2)
mse_df = all_points_df %>% group_by(param,model,nsample,case) %>% 
      summarise(MSE = mean(sq_error), sd = sd(sq_error)) %>% filter(model!='Stan')
all_points_df$nsample = as.factor(all_points_df$nsample)
mse_df$Model = as.factor(mse_df$model)


mse_plot = ggplot(mse_df,aes(x = case, y=MSE,color = Model,shape=Model,
                             linetype = Model)) +  
  geom_point(position = position_dodge(width = 0.5)) + theme_bw() +
  facet_grid( nsample~param,scales = 'free' ) +  
  xlab('Correlation scenario') + ylab('MSE') + 
   geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), 
                 position = position_dodge(width = 0.5),
                 width=.2) +  
  scale_color_manual(values = c("brown1","deepskyblue")) + 
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position = 'bottom')

ggsave(mse_plot,file = paste0('./MSE_plot_',scenario,'.png'),
       width = 6, height = 7)

# Effective sample size per second
sel_points = all_points_df %>% filter(model!='mle') %>% 
  select(model,param,EffSize,nsample,case)
sel_points$rep = rep(rep(rep(1:504,6),2),2)
sel_times = all_times_df %>% filter(model!='mle' )
sel_times$nsample = as.factor(sel_times$nsample)
sel_times$rep = rep(rep(rep(1:504,each = 2),2),2)

eff_df = left_join(sel_points,sel_times)
eff_df = eff_df %>% mutate(eff_per_sec = EffSize/CompTimes)

ess_plot = ggplot(eff_df) + geom_boxplot(aes(x =case, y =eff_per_sec,col=model,
                                             linetype = model)) +
  facet_grid(nsample~param) +theme_bw() + 
  xlab('Correlation scenario') + ylab('Ess/second') + 
  scale_color_manual(values = c("brown1","chartreuse3")) +
  scale_linetype_manual(values = c("Gibbs" = 1, 
                                   "Stan" = 6))+ 
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position = 'bottom')


ggsave(ess_plot,file = paste0('./ESS_plot_',scenario,'.png'),
       width = 6, height = 7)

