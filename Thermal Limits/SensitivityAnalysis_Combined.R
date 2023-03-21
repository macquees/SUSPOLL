## Last edited 2023/03/14 by Sarah MacQueen

# Take the data generated in SensitivityAnalysis3.m and do the regression 
# as described in Saltelli2006

setwd("~/SUSPOLL/Code/Paper 1 code and data to publish")
# install.packages("dismo")
# install.packages("gbm")
library(dismo)
library(gbm)

interaction_level = 2   #level of interactions to consider in the BRT fitting (tree complexity)
VaryEnvironment = TRUE  #always true - vary environmental parameters
Bumblebee = FALSE
CoolingOn = FALSE 

### Read in the data for parameter samples and heat flux ###

if(Bumblebee == TRUE){  #Bumblebee
  if(CoolingOn==TRUE){
    Equilibria = read.csv(file="Thorax_Equilibria_Variability_combined_10000_BB_CoolingOn.csv",header=FALSE)
  }else{
    Equilibria = read.csv(file="Thorax_Equilibria_Variability_combined_10000_BB_CoolingOff.csv",header=FALSE)
  }
  Parameter_Values = read.csv(file="ParameterSample_10000_combined_BB.csv")
  Parameter_Values = Parameter_Values[,-1]      #remove 1st column of indices
  colnames(Parameter_Values) = c("delta_T_h", "i0", "i0_flying", "M_b", "M_th", "E", "c", "r", "T_mK", "alpha_si", 
                                   "epsilon", "A_th", "A_h", "alpha_so", "alpha_th", "f", "P", "T_g", "T_air", "epsilon_e",
                                   "s", "C_l", "n", "l_th", "v", "T_0", "R_0", "D_A", "h_fg", "rh")
    
  if(VaryEnvironment == TRUE){
    Parameter_Values = Parameter_Values[,c(-3,-27,-28,-29,-30)]  
    #remove flying i0, HB only values
    n_params = dim(Parameter_Values)[2] 
  }else{
    Parameter_Values = Parameter_Values[,c(-3,-16,-17,-18,-19,-27,-28,-29,-30)]  
    #remove flying i0, weather values, HB only values
    # colnames(Parameter_Values) = c('delta_T_h','i0','M_b','E','M_th',
    #                                'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
    #                                'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
    #                                'c_l','n','l_th','T_0')
    n_params = dim(Parameter_Values)[2] 
  }
}else{   #Honeybee
  if(CoolingOn==TRUE){
    Equilibria = read.csv(file="Thorax_Equilibria_Variability_combined_30000_HB_CoolingOn.csv",header=TRUE)
    Equilibria = Equilibria[,-1]  #remove 1st column with row indices
    Parameter_Values = read.csv(file="ParameterSample_30000_combined_HB.csv")   #the combined file
    }else{  #cooling off
      Equilibria = read.csv(file="Thorax_Equilibria_Variability_combined_30000_HB_CoolingOff.csv",header=TRUE)
      Equilibria = Equilibria[,-1]  #remove 1st column with row indices
      Parameter_Values = read.csv(file="ParameterSample_30000_combined_HB.csv")  #the original file
    }
  Parameter_Values = Parameter_Values[,-1]      #remove 1st column of row indices
  colnames(Parameter_Values) = c("delta_T_h", "i0", "i0_flying", "M_b", "M_th", "E", "c", "r", "T_mK", "alpha_si", 
                                 "epsilon", "A_th", "A_h", "alpha_so", "alpha_th", "f", "P", "T_g", "T_air", "epsilon_e",
                                 "s", "C_l", "n", "l_th", "v", "T_0", "R_0", "D_A", "h_fg", "rh")
  if(VaryEnvironment == TRUE){
    Parameter_Values = Parameter_Values[,c(-3,-8)]  
    #remove flying i0, BB only values
    n_params = dim(Parameter_Values)[2] 
  }else{
    Parameter_Values = Parameter_Values[,c(-3,-8,-16,-17,-18,-19,-30)]  
    #remove flying i0, weather values, BB only values
    n_params = dim(Parameter_Values)[2] 
  }
}





reg_data_equilibria = cbind(Equilibria,Parameter_Values)
# colnames(reg_data_equilibria) = c('equilibria','delta_T_h','i0','M_b','E','M_th',
#                                           'c','r','T_mK','alpha_si','epsilon_a','A_th','A_h',
#                                           'alpha_s0','alpha_th','a','P','T_gC','T_aC','epsilon_e',
#                                           'c_l','n','l_th','T_0')
reg_data_equilibria_edited = na.omit(reg_data_equilibria)

hist(reg_data_equilibria_edited[ ,1], main='equilibrium thorax temps',xlab='temperature')

# maxouts=which(reg_data_equilibria_edited[ ,1]>=100)
maxouts=which(reg_data_equilibria_edited[ ,1]>=500)

if(length(maxouts)>0){
  reg_data_equilibria_edited = reg_data_equilibria_edited[-maxouts, ]  #remove the rows with temp 100
}
hist(reg_data_equilibria_edited[ ,1], main='equilibrium thorax temps',xlab='temperature')




######################## Fitting #################################

BRT_equilibria = gbm.step(data=reg_data_equilibria_edited,gbm.x=2:(n_params+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
#par(mfrow=c(4,4))
#gbm.plot(BRT_equilibria_flying,n.plots=16,write.title=FALSE)    #plots
contributions = BRT_equilibria$contributions      #the influence of each parameter
R2 = BRT_equilibria$cv.statistics$correlation.mean^2 #might be the R^2 value... 


contributions
R2

### interactions between parameters 
# Interactions_equilibria = gbm.interactions(BRT_equilibria)
# interactions = Interactions_equilibria$interactions
# InteractionsRank = Interactions_equilibria$rank.list

if(Bumblebee==TRUE){
  if(CoolingOn==TRUE){
    write.csv(contributions,file="combined_contributions_bumblebee_coolingon.csv")
    write.csv(R2,file="combined_R2_value_bumblebee_coolingon.csv")
  }else{ #no cooling
    write.csv(contributions,file="combined_contributions_bumblebee_coolingoff.csv")
    write.csv(R2,file="combined_R2_value_bumblebee_coolingoff.csv")
  }
}else{ #honeybee
  if(CoolingOn==TRUE){
    write.csv(contributions,file="combined_contributions_honeybee_coolingon.csv")
    write.csv(R2,file="combined_R2_value_honeybee_coolingon.csv")
  }else{ #no cooling
    write.csv(contributions,file="combined_contributions_honeybee_coolingoff.csv")
    write.csv(R2,file="combined_R2_value_honeybee_coolingoff.csv")
  }
}








######################################################################################################
########### Check for sufficient sample size ###################################

### Check the sampling efficiency
length = dim(reg_data_equilibria_edited)[1]

### sub-samples
#make them a fraction of the total size in case of reduction due to maxouts
sub_10 = floor(length*.1)
sub_25 = floor(length*.25)
sub_50 = floor(length*.5)
sub_75 = floor(length*.75)

sample_1000 = sample(1:length,size=sub_10,replace=FALSE)
sample_2500 = sample(1:length,size=sub_25,replace=FALSE)
sample_5000 = sample(1:length,size=sub_50,replace=FALSE)
sample_7500 = sample(1:length,size=sub_75,replace=FALSE)


### BRTs for each sub-sample
reg_data_equilibria_edited_1000 = reg_data_equilibria_edited[sample_1000,]
BRT_1000 = gbm.step(data=reg_data_equilibria_edited_1000,gbm.x=2:(n_params+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_1000 = BRT_1000$contributions[,2]      #the influence of each parameter

reg_data_equilibria_edited_2500 = reg_data_equilibria_edited[sample_2500,]
BRT_2500 = gbm.step(data=reg_data_equilibria_edited_2500,gbm.x=2:(n_params+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_2500 = BRT_2500$contributions[,2]      #the influence of each parameter

reg_data_equilibria_edited_5000 = reg_data_equilibria_edited[sample_5000,]
BRT_5000 = gbm.step(data=reg_data_equilibria_edited_5000,gbm.x=2:(n_params+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_5000 = BRT_5000$contributions[,2]      #the influence of each parameter

reg_data_equilibria_edited_7500 = reg_data_equilibria_edited[sample_7500,]
BRT_7500 = gbm.step(data=reg_data_equilibria_edited_7500,gbm.x=2:(n_params+1),gbm.y=1,learning.rate=0.01, bag.fraction=0.75, tree.complexity = interaction_level, n.folds=10, family="gaussian")
Influences_7500 = BRT_7500$contributions[,2]      #the influence of each parameter

Influences_10000 = BRT_equilibria$contributions[,2]


### calculate the index D 
# For 1000-2500
Average_Influences_1000_2500 = apply(cbind(Influences_1000,Influences_2500),1,mean)
Inner_sum_1000 = sum(Influences_1000*log(Influences_1000)/2,na.rm=TRUE)
Inner_sum_2500 = sum(Influences_2500*log(Influences_2500)/2,na.rm=TRUE)
D_1000_2500 = Inner_sum_1000+Inner_sum_2500 - sum(Average_Influences_1000_2500*log(Average_Influences_1000_2500),na.rm=TRUE)


# For 2500-5000
Average_Influences_2500_5000 = apply(cbind(Influences_2500,Influences_5000),1,mean)
Inner_sum_2500 = sum(Influences_2500*log(Influences_2500)/2,na.rm=TRUE)
Inner_sum_5000 = sum(Influences_5000*log(Influences_5000)/2,na.rm=TRUE)
D_2500_5000 = Inner_sum_2500+Inner_sum_5000 - sum(Average_Influences_2500_5000*log(Average_Influences_2500_5000),na.rm=TRUE)


# For 5000-7500 
Average_Influences_5000_7500 = apply(cbind(Influences_5000,Influences_7500),1,mean)
Inner_sum_5000 = sum(Influences_5000*log(Influences_5000)/2,na.rm=TRUE)
Inner_sum_7500 = sum(Influences_7500*log(Influences_7500)/2,na.rm=TRUE)
D_5000_7500 = Inner_sum_5000+Inner_sum_7500 - sum(Average_Influences_5000_7500*log(Average_Influences_5000_7500),na.rm=TRUE)


# For 7500-10000
Average_Influences_7500_10000 = apply(cbind(Influences_7500,Influences_10000),1,mean)
Inner_sum_7500 = sum(Influences_7500*log(Influences_7500)/2,na.rm=TRUE)
Inner_sum_10000 = sum(Influences_10000*log(Influences_10000)/2,na.rm=TRUE)
D_7500_10000 = Inner_sum_7500+Inner_sum_10000 - sum(Average_Influences_7500_10000*log(Average_Influences_7500_10000),na.rm=TRUE)



### All together
D = exp(c(D_1000_2500,D_2500_5000,D_5000_7500,D_7500_10000))

par(mfrow=c(1,1))
plot(D,main="Honeybee, Cooling Off",type='b',
     ylab = "D_beta")


### save the results
if(Bumblebee==TRUE){
  if(CoolingOn==TRUE){
    write.csv(Influences_1000,file="combined_Influences_1000_bumblebee_coolingon.csv")
    write.csv(Influences_2500,file="combined_Influences_2500_bumblebee_coolingon.csv")
    write.csv(Influences_5000,file="combined_Influences_5000_bumblebee_coolingon.csv")
    write.csv(Influences_7500,file="combined_Influences_7500_bumblebee_coolingon.csv")
    write.csv(Influences_10000,file="combined_Influences_10000_bumblebee_coolingon.csv")
    write.csv(D,file="combined_D_values_bumblebee_coolingon.csv")
  }else{ #no cooling
    write.csv(Influences_1000,file="combined_Influences_1000_bumblebee_coolingoff.csv")
    write.csv(Influences_2500,file="combined_Influences_2500_bumblebee_coolingoff.csv")
    write.csv(Influences_5000,file="combined_Influences_5000_bumblebee_coolingoff.csv")
    write.csv(Influences_7500,file="combined_Influences_7500_bumblebee_coolingoff.csv")
    write.csv(Influences_10000,file="combined_Influences_10000_bumblebee_coolingoff.csv")
    write.csv(D,file="combined_D_values_bumblebee_coolingoff.csv")
  }
}else{ #honeybee
  if(CoolingOn==TRUE){
    write.csv(Influences_1000,file="combined_Influences_1000_honeybee_coolingon.csv")
    write.csv(Influences_2500,file="combined_Influences_2500_honeybee_coolingon.csv")
    write.csv(Influences_5000,file="combined_Influences_5000_honeybee_coolingon.csv")
    write.csv(Influences_7500,file="combined_Influences_7500_honeybee_coolingon.csv")
    write.csv(Influences_10000,file="combined_Influences_10000_honeybee_coolingon.csv")
    write.csv(D,file="combined_D_values_honeybee_coolingon.csv")
  }else{ #no cooling
    write.csv(Influences_1000,file="combined_Influences_1000_honeybee_coolingoff.csv")
    write.csv(Influences_2500,file="combined_Influences_2500_honeybee_coolingoff.csv")
    write.csv(Influences_5000,file="combined_Influences_5000_honeybee_coolingoff.csv")
    write.csv(Influences_7500,file="combined_Influences_7500_honeybee_coolingoff.csv")
    write.csv(Influences_10000,file="combined_Influences_10000_honeybee_coolingoff.csv")
    write.csv(D,file="combined_D_values_honeybee_coolingoff.csv")
  }
}

