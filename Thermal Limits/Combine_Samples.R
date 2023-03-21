## Last edited 2023/03/14 by Sarah MacQueen
## used to combine the parameter samples and correspOnding equilibria for the honeybee/cooling On and Off cases
## because the original sample size wasn't large enough 

setwd("~/SUSPOLL/Code/Paper 1 code and data to publish")

Parameter_Values_1 = read.csv(file="ParameterSample_10000_combined_HB_1.csv")
Parameter_Values_1 = Parameter_Values_1[,-1]      #remove 1st column of indices
colnames(Parameter_Values_1) = c("delta_T_h", "i0", "i0_flying", "M_b", "M_th", "E", "c", "r", "T_mK", "alpha_si", 
                               "epsilon", "A_th", "A_h", "alpha_so", "alpha_th", "f", "P", "T_g", "T_air", "epsilon_e",
                               "s", "C_l", "n", "l_th", "v", "T_0", "R_0", "D_A", "h_fg", "rh")

Parameter_Values_2 = read.csv(file="ParameterSample_10000_combined_HB_2.csv")
Parameter_Values_2 = Parameter_Values_2[,-1]      #remove 1st column of indices
colnames(Parameter_Values_2) = c("delta_T_h", "i0", "i0_flying", "M_b", "M_th", "E", "c", "r", "T_mK", "alpha_si", 
                                 "epsilon", "A_th", "A_h", "alpha_so", "alpha_th", "f", "P", "T_g", "T_air", "epsilon_e",
                                 "s", "C_l", "n", "l_th", "v", "T_0", "R_0", "D_A", "h_fg", "rh")

Parameter_Values_3 = read.csv(file="ParameterSample_10000_combined_HB_3.csv")
Parameter_Values_3 = Parameter_Values_3[,-1]      #remove 1st column of indices
colnames(Parameter_Values_3) = c("delta_T_h", "i0", "i0_flying", "M_b", "M_th", "E", "c", "r", "T_mK", "alpha_si", 
                                 "epsilon", "A_th", "A_h", "alpha_so", "alpha_th", "f", "P", "T_g", "T_air", "epsilon_e",
                                 "s", "C_l", "n", "l_th", "v", "T_0", "R_0", "D_A", "h_fg", "rh")

Parameter_Values_30000 = rbind(Parameter_Values_1,Parameter_Values_2,Parameter_Values_3)
write.csv(Parameter_Values_30000,file="ParameterSample_30000_combined_HB.csv")

# Parameter_Values_20000 = rbind(Parameter_Values_1,Parameter_Values_2)
# write.csv(Parameter_Values_20000,file="ParameterSample_20000_combined_HB.csv")


### Honeybee, Cooling On
Equilibria_CoolingOn_1 = read.csv(file="Thorax_Equilibria_Variability_combined_10000_HB_CoolingOn_1.csv",header=FALSE)
Equilibria_CoolingOn_2 = read.csv(file="Thorax_Equilibria_Variability_combined_10000_HB_CoolingOn_2.csv",header=FALSE)
Equilibria_CoolingOn_3 = read.csv(file="Thorax_Equilibria_Variability_combined_10000_HB_CoolingOn_3.csv",header=FALSE)

Equilibria_CoolingOn_30000 = rbind(Equilibria_CoolingOn_1,Equilibria_CoolingOn_2,Equilibria_CoolingOn_3)
colnames(Equilibria_CoolingOn_30000) = c("thorax_temp")

write.csv(Equilibria_CoolingOn_30000,file="Thorax_Equilibria_Variability_combined_30000_HB_CoolingOn.csv")


### Honeybee, Cooling Off
Equilibria_CoolingOff_1 = read.csv(file="Thorax_Equilibria_Variability_combined_10000_HB_CoolingOff_1.csv",header=FALSE)
Equilibria_CoolingOff_2 = read.csv(file="Thorax_Equilibria_Variability_combined_10000_HB_CoolingOff_2.csv",header=FALSE)
Equilibria_CoolingOff_3 = read.csv(file="Thorax_Equilibria_Variability_combined_10000_HB_CoolingOff_3.csv",header=FALSE)

Equilibria_CoolingOff_30000 = rbind(Equilibria_CoolingOff_1,Equilibria_CoolingOff_2,Equilibria_CoolingOff_3)
colnames(Equilibria_CoolingOff_30000) = c("thorax_temp")

write.csv(Equilibria_CoolingOff_30000,file="Thorax_Equilibria_Variability_combined_30000_HB_CoolingOff.csv")
