% Last edited 2023/03/14 by Sarah MacQueen
% does SA model runs for a single combined metabolic (and flight speed) range instead
% of resting/flying
% Note: ode15s is required for Honeybee/CoolingOn, but ode45 is fine for
% the other cases

%% 
%%%%%%%%%%% Model Switches %%%%%%%%%%%%%%%%
plotting = false;   %should we plot the solution curves


Bumblebee = false;    %only one of Bumblebee and Honeybee can be true
Honeybee = true;

CoolingOff = true;   %only one of cooling can be true
CoolingSwitch = false;  %not used
CoolingOn = false; 

%unused
% Resting = false;      %only one of metabolic states can be true 
% Shivering = false; 
% Flying = true;
% indicator = 3;  %is the bee resting/shivering/flying

VaryEnvironment = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Code Structure Stuff %%%%%%%%
n_samples = 10000;
cutoffTemp = 500;  %upper limit temp to cut numerical solving at
maxy=cutoffTemp+273.15;
pcolors=['b' 'y' 'r'];
tspan = 0:6000;  %the time over which I have data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Fixed Values for HB and BB %%%%%%%%%%%%%%%%%%%%%%%
Pr = 1.013*10^5; %atmospheric pressure in N/m^2 (value used in Sidebotham)
R_specific = 287.058;  %J/kg/K for dry air
A = 9.1496*10^10;  %clausius-clapyron constant A for water in N/m^2
B = -5.1152*10^3;  %clausius-clapyron constant B for water in K
MM_air = 0.0289652;   %molar mass of dry air in g/mol
MM_vapor = 0.018016;  %molar mass of water vapor in kg/mol
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant
k = 8.617333262145*10^(-5);   %Bolzmann's constant


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Storage Vectors %%%%%%%%%%%%
Thorax_Equilibria_Variability = zeros(n_samples,1);

if Honeybee == true
    Parameter_Values = readtable('ParameterSample_10000_combined_HB_3.csv','ReadVariableNames',true);
    %row = sample; column = parameter (1st column is row numbers...)
    %_1 is original sample, _2 and _3 are additional samples 
    %no number is combination of _1, _2, _3
    masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
end

if Bumblebee == true
    Parameter_Values = readtable('ParameterSample_10000_combined_BB.csv','ReadVariableNames',true);
    %row = sample; column = parameter (1st column is row numbers...)
    masses = [0.177 0.177 0.177];   %reference weight for Kammer (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
end

%prepare the plot for ploting solution curves and histogram 
tplot= tiledlayout(1,1);   %Create a 1-by-1 tiled chart layout tplot.
ax1 = axes(tplot);   %Create an axes object ax1 by calling the axes function and specifying tplot as the parent object.

for sample_index=1:n_samples
% for sample_index=9999
sample_index  %print to show where we are

%%%%% Bee Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_T_h = Parameter_Values.V1(sample_index);
i_0 = Parameter_Values.V2(sample_index);
% I_resting = Parameter_Values.V2(sample_index);  %note to check T_ref if using Heinrich vs Kammer data
% I_flying = Parameter_Values.V3(sample_index);  %note to check T_ref if using Heinrich vs Kammer data
M_b = Parameter_Values.V4(sample_index);   %mass of the bee in g (the default used by Cooper1985 is 100mg)
E = Parameter_Values.V5(sample_index);    %Brown2004 activation energy
M_th = Parameter_Values.V6(sample_index);  %mass of thorax in g (mean reported in Cooper1985)
c = Parameter_Values.V7(sample_index);  %specific heat (in cal/g*degC converted to J/g*degC), cited in May1976
r = Parameter_Values.V8(sample_index);   %abdomen cooling parameter
T_mK = Parameter_Values.V9(sample_index);  %midpoint of cooling switch function
alpha_si = Parameter_Values.V10(sample_index);     %shape factor for incoming solar radiation (Cooper1985)
epsilon_a = Parameter_Values.V11(sample_index);   % absorptivity of bees (Willmer1981)
A_th = Parameter_Values.V12(sample_index)  ; %in m^2, from Cooper1985
A_h = Parameter_Values.V13(sample_index)  ; %head surface area in m^2, from Cooper1985
alpha_so = Parameter_Values.V14(sample_index);     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = Parameter_Values.V15(sample_index);     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
epsilon_e = Parameter_Values.V20(sample_index);       %(fill in the reference for this!)
s = Parameter_Values.V21(sample_index);   %fraction of core temp at surface
C_l = Parameter_Values.V22(sample_index);   %(fill in reference) 
n = Parameter_Values.V23(sample_index);       %(fill in reference)
l_th = Parameter_Values.V24(sample_index);   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
v = Parameter_Values.V25(sample_index); %flight speed or wind speed (3.1m/s from Cooper1985)
y0 = Parameter_Values.V26(sample_index)+273.15;   %initial temperature of the bee's head in K (observed by me! ish)
R_0 = Parameter_Values.V27(sample_index);  %radius of nectar droplet, 1mm in m
D_A = Parameter_Values.V28(sample_index);  %diffusion coefficient of air into itself in m^2/s -
h_fg = Parameter_Values.V29(sample_index);  %latent heat of vaporization of water, J/kg 

if VaryEnvironment == true
    P = Parameter_Values.V17(sample_index); %mean solar irradiance from all of Arrian's data
    T_gC = Parameter_Values.V18(sample_index);                  %ground surface temp in C https://www.met.ie/climate/available-data/monthly-data %Phoenix Park June 2021
    T_aC = Parameter_Values.V19(sample_index);    %mean air temp in C from Arrian's data
    rh =  Parameter_Values.V30(sample_index);   %mean relative humidity from Arrian's data
    a = Parameter_Values.V16(sample_index);   %called f in paper, fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)
else
    P = 332.3878; %mean solar irradiance from all of Arrian's data
    T_gC = 17.1;                  %ground surface temp in C https://www.met.ie/climate/available-data/monthly-data %Phoenix Park June 2021
    T_aC = 20;    %air temp
    rh =  0.6908;   %mean relative humidity from Arrian's data
    a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)
end
T_aK = T_aC+273.15;     %air temp in K
T_gK = T_gC+273.15;        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report


%%%%%  Environment Parameters that depend on air temp %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = (0.02646*T_aK^1.5)/(T_aK+245.4*10^(-12/T_aK));
mu = (1.458*10^(-6)*T_aK^1.5)/(T_aK+110.4); %dynamic viscosity of air 
rho = Pr/(R_specific*T_aK);  %in humid air air at sea level 
nu = mu/rho; %kinematic viscosity %The Shock Absorber Handbook



%%%%%% Solar Radiation %%%%%%%
% Does not depend on T_th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = (alpha_si*epsilon_a*A_th*P)/(M_th*c) + (alpha_so*epsilon_a*A_th*a*P)/(M_th*c);   %solar radiation thorax
S_h = (alpha_si*epsilon_a*A_h*P)/(M_th*c) + (alpha_so*epsilon_a*A_h*a*P)/(M_th*c);   %solar radiation head

%%%%%% Thermal Radiation %%%%%%%
% Term 2 does depend on T_th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1 = (alpha_th*epsilon_a*A_th*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee thorax (from sky + earth)
R2 = (epsilon_e*A_th*sigma)/(M_th*c);   %thermal radiation out of bee thorax
R1_h = (alpha_th*epsilon_a*A_h*(delta*T_aK.^6+sigma*T_gK.^4))/(M_th*c);  %thermal radiation into bee head (from sky + earth)
R2_h = (epsilon_e*A_h*sigma)/(M_th*c);   %thermal radiation out of bee head


%%%%%%%% Convection %%%%%%%%%
%1st term does depend on T_th
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = (C_l*kappa/l_th)*(v*l_th/nu)^n;
C1 = (h*A_th*s)/(M_th*c);           %thorax, will be multiplied by T_th, *-.9 is for thorax surface temperature
C2 = (-h*A_th*T_aK)/(M_th*c);           %thorax
C1_h = (h*A_h)/(M_th*c);           %head
C2_h = (-h*A_h*T_aK)/(M_th*c);           %head

%%%%%%%%%%% Metabolic %%%%%%%%%
%Does  depend on T_th 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_ref = masses(1);
T_ref = RefTemps(1);
% norm_constants = [I_resting, I_flying, I_flying];   %resting/thermoregulating/flying = 1,2,3
% i_0 = norm_constants(indicator);

I = (i_0*((M_b/M_ref)^(3/4))*exp((E/(k*T_ref))))*(1/(M_th*c)); %uses i0=I, depends on T_th, relative to ref temp/mass



%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatflux.m, which may need to be
%renamed

%the times at which the data for the time dependent funcions was measured
%%%%%%%%%%%%% Solve ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note: the ODE function is defind in heatflux.m, which may need to be
%renamed
Opt = odeset('Events',@(t,y)myEvent(t,y,maxy),'RelTol',1e-12,'AbsTol',1e-14);   %using this will make ode45 stop solving at T_th=maxy






%for now, only set for using cooling switching 
Cooling_vals = Cooling_Flux(Bumblebee,r,M_th,c,A,B,T_aK,rh,Pr,MM_air,MM_vapor,R_specific,R_0,D_A,h_fg);
Ab = Cooling_vals(1);
Ev1 = Cooling_vals(2);
Ev2 = Cooling_vals(3);
CoolingSwitch_indicator = 1;
[t,y] = ode15s(@(t,y) heatfluxhead(t,y,S,R1,R2,C1,C2,I,S_h,R1_h,R2_h,C1_h,C2_h,delta_T_h,E,k,T_aK,T_mK,A,B,rh,Pr,MM_air,MM_vapor,Honeybee,CoolingSwitch_indicator,Ab,Ev1,Ev2), tspan, y0,Opt); %for use with constant environmental conditions
if plotting ==true
plot(t,y-273.15,'color',[0,0,0]+0.5);   %plot in celsius  'r' makes a red line; 'b' = blue; 'k' = black
hold on;
end
y = calc_length(tspan,y,maxy);  %check if bee maxed out temp or got cut off with imaginary parts and fill in if necessary

if max(y) == cutoffTemp+273.15   %if temp reached cutoff temp, use that
    Thorax_Equilibria_Variability(sample_index) = cutoffTemp;
else  %otherwise, take the mean/median of the end
    Thorax_Equilibria_Variability(sample_index) = median(y(5500:length(tspan)))-273.15;     %median instead of mean should be more stable against cycles
end

    
end

if Honeybee==true
%     if indicator==3
%         writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_flying_combined_10000_HB_CoolingOn.csv');
%     end
%     if indicator==2
%         writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_shivering_combined_10000_HB_CoolingOn.csv');
%     end
%     if indicator==1
%         writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_resting_combined_10000_HB_CoolingOn.csv');
%     end
    if CoolingOn == true
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_combined_10000_HB_CoolingOn.csv');
    else
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_combined_10000_HB_CoolingOff_3.csv');
    end

end

if Bumblebee==true
%     if indicator==3
%         writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_flying_combined_10000_BB_CoolingOn.csv');
%     end
%     if indicator==2
%         writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_shivering_combined_10000_BB_CoolingOn.csv');
%     end
%     if indicator==1
%         writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_resting_combined_10000_BB_CoolingOn.csv');
%     end
    if CoolingOn == true
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_combined_10000_BB_CoolingOn.csv');
    else
        writematrix(Thorax_Equilibria_Variability,'Thorax_Equilibria_Variability_combined_10000_BB_CoolingOff.csv');
    end
end

if plotting ==true
%add the labels to the solution curve plots
title('Flying Bee');
subtitle(['P=',num2str(P),'W/m^2'])
xlabel('Time (s)') ;
ylabel('Thorax Temperature (C)') ;
if Bumblebee == true
yline(30,'color','r');  %30 is the min thorax temp for flight, Heinrich1983
yline(45,'color','r');   %lethal thorax temp Heinrich1976?ish
end
if Honeybee == true
yline(35,'color','r');  %30 is the min thorax temp for flight, Heinrich1983
yline(52,'color','r');   %lethal thorax temp Heinrich1976?ish
end
ylim([0,cutoffTemp])    %set the y-axis limits to match the rotated histogram
end


% %%%%% Plot a nice sideways histogram on top of it! %%%%%%
% ax2 = axes(tplot);  %Create an axes object ax2 by calling the axes function and specifying t as the parent object.
% h = histogram(Thorax_Equilibria_Variability,'orientation','horizontal');  %establish the histogram plot
% set(get(h,'Parent'),'xdir','r')   %put it on the other axis
% ax2.XAxisLocation = 'top'; %Move the x-axis to the top,  
% ax2.YAxisLocation = 'right';   %and move the y-axis to the right.
% ylim([0,cutoffTemp])    %set the y-axis limits for the rotated histogram
% set(ax2(1),'XTickLabel','');   %but don't actually show the axis labels
% set(ax2(1),'YTickLabel','');
% xlim([0,2000])  %reset the x-limits (for the histogram only), for visual effect
% ax2.Color = 'none';  %Set the color of the axes object to 'none' so that the underlying plot is visible.
% ax1.Box = 'off';    %Turn off the plot boxes to prevent the box edges from obscuring the x- and y-axes. 
% ax2.Box = 'off';
% 

hold off 



beep
