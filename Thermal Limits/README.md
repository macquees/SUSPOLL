Code to analyze the model and produce all results and figures in the paper Thermal limits of bumblebees and honeybees are modulated by different functional traits: predictions of a mechanistic model.
Equilibrium thorax temperature as a function of air temperature

Used to produce figure 2

The code has to be run separately for bumblebees and honeybees. This is controlled at the top of the file LookUpTable_Script.m, in the “Model Switches” section. 

    • LookupTable_Script.m  produces the data used to make figure 2 by running the model for various behaviour combinations and a range of air temperatures from 0 to 50 C. 
        ◦ Calls: Parameters_Script.m which contains the default parameter values 
        ◦ Uses: myEvent.m function, which is used to stop integration at a specified thorax temperature (upper limit)
        ◦ Calls: SolveModel_Script.m which uses a built in ode solver to solve the model
            ▪ Uses: heatfluxhead.m function, which contains the differential equation for the model 
            ▪ Uses: Cooling_Flux.m function, which calculates the abdomen or nectar cooling values
        ◦ Uses: calc_length.m function, which fills in t and y values for the full length of tspan if the solution maxed out
        ◦ Outputs: Flying_Honeybee_LookupTable_Full.csv and Flying_Bumblebee_LookupTable_Full.csv containing the full results as used in the supporting information
        ◦ Outputs: Flying_Honeybee_LookupTable.csv and Flying_Bumblebee_LookupTable.csv containing the full results the results as used in figure 2 of the paper
    • SmoothPlots.R is used to produce the plots for the figures, and the accompanying legends where relevant. 
        ◦ Reads in: Flying_Bumblebee_LookupTable_Full.csv and Flying_Honeybee_LookupTable_Full.csv containing the full results as used in the supporting information
        ◦ Reads in: Flying_Bumblebee_LookupTable.csv, Flying_Honeybee_LookupTable.csv, Resting_Bumblebee_LookupTable.csv,  and Resting_Honeybee_LookupTable.csv  containing the results as used in figure 2 of the paper

Maximum thorax temperature for sustained flight as a function of thorax mass and flight speed

Used to produce figure 3  

The code has to be run separately for bumblebees and honeybees. This is controlled at the top of the file MakePlots.m, in the “Setup” section. The file is run only for a flying bee with cooling behaviour (CoolingOn= true and Flying = true), but the other options (set to false) are included so that the code works.
    • Mass_FlightSpeed.m runs the model for a range of different thorax masses and three different flight speeds. For each thorax mass/flight speed combination, the equilibrium thorax is calculated for air temperatures from 0 to 50 C and the highest air temperature at which flight can be sustained (the bee doesn’t overheat) is found. Both bumblebees and honeybees are included automatically.
        ◦ Calls Parameters_Script.m to load the default parameter values
        ◦ Calls: SolveModel_Script.m which uses a built in ode solver to solve the model
            ▪ Uses: heatfluxhead.m function, which contains the differential equation for the model 
            ▪ Uses: Cooling_Flux.m function, which calculates the abdomen or nectar cooling values
        ◦ Uses: calc_length.m function, which fills in t and y values for the full length of tspan if the solution maxed out
        ◦ Outputs MaxAirTemp_BB.csv and MaxAirTemp_HB.csv which contain the maximum air temperature for sustained flight with a row for each mass and a column for each flight speed. 
    • SmoothPlots.R is used to produce the plots for the figure 
        ◦ Reads in: MaxAirTemp_BB.csv and MaxAirTemp_HB.csv containing the full results the results as used in figure 2 of the paper

Fitting i_0 and E

    • GridSearchFit_General.m evaluates the model at each point on a grid of i_0 and E values. It solves for the air temperature at which the bee is warm enough to start flying and the air temperature at which the bee has to start using its behavioural cooling mechanism. The file runs first bumblebees, then honeybees.
        ◦ Uses: RunModelABC_v2.m function, which takes the grid values for i_0 and E and runs the model 
            ▪ Uses: myEvent.m function, which is used to stop integration at a specified thorax temperature (upper limit)
            ▪ Uses: heatfluxhead.m function, which contains the differential equation for the model 
            ▪ Uses: calc_length.m function, which fills in t and y values for the full length of tspan if the solution maxed out
        ◦ Outputs: Fly_P1_P2_Bumblebee.csv, Cool_P1_P2_Bumblebee.csv, Fly_P1_P2_Honeybee.csv, Cool_P1_P2_Honeybee.csv, which contain a grid of air temperature values corresponding to i_0 (row) and E (column) values
    • GridSearchInterpolation.R runs the simplified Approximate Bayesian Computation fitting procedure. Only one of bumblebee/honeybee can be run at a time. This is set at the top of the file. 
        ◦ Reads in: Fly_P1_P2_Bumblebee.csv, Cool_P1_P2_Bumblebee.csv, Fly_P1_P2_Honeybee.csv, Cool_P1_P2_Honeybee.csv
        ◦ Produces: heatmaps and histogram figure used in Supporting Information
        ◦ Produces files containing the data for the heatmaps and histograms: 
            ▪ Fit_CoolTemp_values_list_BB.csv
            ▪ Fit_CoolTemp_values_list_HB.csv
            ▪ Fit_Dist_values_list_BB.csv
            ▪ Fit_Dist_values_list_HB.csv
            ▪ Fit_E_values_list_BB.csv
            ▪ Fit_E_values_list_HB.csv
            ▪ Fit_FlyTemp_values_list_BB.csv
            ▪ Fit_FlyTemp_values_list_HB.csv
            ▪ Fit_i0_values_list_BB.csv
            ▪ Fit_i0_values_list_HB.csv
        ◦ The i_0 and E fitted values are i0_guess_minvalue and E_guess_minvalue. These may vary slightly with each run of the file due to randomness in the procedure. 


Sensitivity Analysis 

The sensitivity analysis has to be run separately for bumblebees and honeybees, and for abdomen/nectar behavioural cooling on or off. These are controlled at the top of the file SensitivityAnalysis3.m, in the “Model Switches” section. 
    • VariabilityDueToParameters.R is used to generate the parameter sample. It must be run once for bumblebees and once for honeybees, controlled with the logicals at the beginning of the file. 
        ◦ Outputs the files ParameterSample_10000_combined_BB.csv and ParameterSample_10000_combined_HB.csv for bumblebee and honeybee respectively. 
    • SensitivityAnalysis3.m runs the model (calculates the equilibrium thorax temperature) on each parameter sample 
        ◦ Reads in:  ParameterSample_10000_combined_BB.csv or ParameterSample_10000_combined_HB.csv 
        ◦ Uses: myEvent.m function, which is used to stop integration at a specified thorax temperature (upper limit)
        ◦ Uses: Cooling_Flux.m function, which calculates the abdomen or nectar cooling values
        ◦ Uses: heatfluxhead.m function, which contains the differential equation for the model 
        ◦ Outputs files of the form “Thorax_Equilibria_Variability_combined_10000_beetype_coolingtype.csv” where bee types are BB and HB, and cooling types are CoolingOn and CoolingOff. These contain the equilibrium thorax temperature for each parameter sample. 
    • Combine_Samples.R combines multiple latin hypercube parameter sample sets and corresponding thorax equilibria for honeybees
        ◦ Reads in  ParameterSample_10000_combined_HB_sampleset.csv for sample sets 1, 2, and/or 3 (each with 10000 samples)
        ◦ Reads in Thorax_Equilibria_Variability_combined_10000_HB_coolingtype_sampleset.csv for sample sets 1, 2, and/or 3 (each with 10000 samples) and cooling on/off
        ◦ Outputs ParameterSample_numberofsamples_combined_HB.csv Thorax_Equilibria_Variability_combined_numberofsamples_HB_coolingtype.csv for the number of samples combined and cooling on/off
    • SensitivityAnalysis_Combined.R runs the sensitivity analysis: calculates the relative contribution of each parameter, the interactions between parameters, and checks the sampling sufficiency 
        ◦ Bee type and cooling type are set at the beginning of the file – must be run for all 4 combinations separately
        ◦ Reads in ParameterSample_numberofsamples_combined_HB.csv  and Thorax_Equilibria_Variability_combined_numberofsamples_beetype_coolingtype.csv files depending on the bee, cooling type, and number of samples
        ◦ Outputs files of the form: 
            ▪ combined_contributions_beetype_coolingtype.csv
            ▪ combined_D_values_beetype_coolingtype.csv
            ▪ combined_Influences_subsample_beetype_coolingtype.csv
            ▪ combined_R2_value_beetype_coolingtype.csv
        ◦ Produces plots

Nondimensionalisation
    • Bee_model_nodims_v3.nb nondimensionalises a simpler version of the model 
        ◦ no heat transfer to the head (a Q_h term is given, but not included in the model for nondimensionalisation)
        ◦ no behavioural cooling mechanisms (abdomen/nectar)

This workhas emanated from research conducted with the financial support of Science Foundation Ireland as part of the SUSPOLL project under Grant number 17/CDA/4689.
The opinions, findings and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the Science Foundation Ireland.
