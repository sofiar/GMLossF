####################################################################
###################### Get STAN results ############################ 
####################################################################

# T = 30
load(file='/u/ruizsuar/GMLossF/Main_scripts/StanNsamp30.RData')

results_b_30 = matrix(NA,ncol = 2, nrow =3)
colnames(results_b_30) = c('triplet_Bay','triplet_freq')
rownames(results_b_30) = c('Bias','Variace','MSE')

results_theta1_30 = matrix(NA,ncol = 2, nrow =3)
colnames(results_theta1_30) = c('triplet_Bay','triplet_freq')
rownames(results_theta1_30) = c('Bias','Variace','MSE')

results_theta2_30 = matrix(NA,ncol = 2, nrow =3)
colnames(results_theta2_30) = c('triplet_Bay','triplet_freq')
rownames(results_theta2_30) = c('Bias','Variace','MSE')

for (s in 1:s){

    # b
    results_b_30[1,s] = mean(output_b_mean[,s]) - b_real[s] # bias
    results_b_30[2,s] = mean((output_b_mean[,s] - mean(output_b_mean[,s]))^2) # variance
    results_b_30[3,s] = mean((output_b_mean[,s] - b_real[s])^2) # MSE

    # theta1
    results_theta1_30[1,s] = mean(output_theta1_mean[,s]) - theta1_real[s] # bias
    results_theta1_30[2,s] = mean((output_theta1_mean[,s] - mean(output_theta1_mean[,s]))^2) # variance
    results_theta1_30[3,s] = mean((output_theta1_mean[,s] - theta1_real[s])^2) # MSE

    # theta2
    results_theta2_30[1,s] = mean(output_theta2_mean[,s]) - theta2_real[s] # bias
    results_theta2_30[2,s] = mean((output_theta2_mean[,s] - mean(output_theta2_mean[,s]))^2) # variance
    results_theta2_30[3,s] = mean((output_theta2_mean[,s] - theta2_real[s])^2) # MSE
}

results_b_30
results_theta1_30
results_theta2_30


# T = 100
load(file='/u/ruizsuar/GMLossF/Main_scripts/StanNsamp100.RData')

results_b_100 = matrix(NA,ncol = 2, nrow =3)
colnames(results_b_100) = c('triplet_Bay','triplet_freq')
rownames(results_b_100) = c('Bias','Variace','MSE')

results_theta1_100 = matrix(NA,ncol = 2, nrow =3)
colnames(results_theta1_100) = c('triplet_Bay','triplet_freq')
rownames(results_theta1_100) = c('Bias','Variace','MSE')

results_theta2_100 = matrix(NA,ncol = 2, nrow =3)
colnames(results_theta2_100) = c('triplet_Bay','triplet_freq')
rownames(results_theta2_100) = c('Bias','Variace','MSE')

for (s in 1:s){

    # b
    results_b_100[1,s] = mean(output_b_mean[,s]) - b_real[s] # bias
    results_b_100[2,s] = mean((output_b_mean[,s] - mean(output_b_mean[,s]))^2) # variance
    results_b_100[3,s] = mean((output_b_mean[,s] - b_real[s])^2) # MSE

    # theta1
    results_theta1_100[1,s] = mean(output_theta1_mean[,s]) - theta1_real[s] # bias
    results_theta1_100[2,s] = mean((output_theta1_mean[,s] - mean(output_theta1_mean[,s]))^2) # variance
    results_theta1_100[3,s] = mean((output_theta1_mean[,s] - theta1_real[s])^2) # MSE

    # theta2
    results_theta2_100[1,s] = mean(output_theta2_mean[,s]) - theta2_real[s] # bias
    results_theta2_100[2,s] = mean((output_theta2_mean[,s] - mean(output_theta2_mean[,s]))^2) # variance
    results_theta2_100[3,s] = mean((output_theta2_mean[,s] - theta2_real[s])^2) # MSE
}

results_b_100
results_theta1_100
results_theta2_100

