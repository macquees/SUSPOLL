%%%%%%%%%%%%%%%%%  Contains all default parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcolors=['b' 'y' 'r' 'g' 'k'];   %plotting colors


%%%%%% Constant values %%%%%%%%
k = 8.617333262145*10^(-5);   %Bolzmann's constant
delta = 5.31*10^(-13);   %fill in the name of this constant
sigma = 5.67*10^(-8);   %fill in the name of this constant
h_fg = 2.3819*10^6;  %latent heat of vaporization of water at 50C, J/kg  (it actually depends on temp)
%https://www.engineeringtoolbox.com/water-properties-d_1573.html
A = 9.1496*10^10;  %clausius-clapyron constant A for water in N/m^2
B = -5.1152*10^3;  %clausius-clapyron constant B for water in K
rh =  0.6908;   %mean relative humidity from Arrian's data
%D_A = 2.06*10^-5;  %diffusion coefficient of air into itself in m^2/s - %calculated version is below
MM_air = 0.0289652;   %molar mass of dry air in kg/mol
%https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
MM_vapor = 0.018016;  %molar mass of water vapor in kg/mol

%%%%%%%%%%% Environmental Parameters %%%%%%%%%%%%%%%%
T_aC = 20;    %air temp
T_aK = T_aC+273.15;        %air surface temp in K
P = 332.3878; %mean solar irradiance from all of Arrian's data
T_gC = 17.1;                  %ground surface temp in C https://www.met.ie/climate/available-data/monthly-data %Phoenix Park June 2021
T_gK = T_gC+273.15;        %ground surface temp in K, very vague estimate from https://www.met.ie/forecasts/farming/agricultural-data-report
Pr = 1.013*10^5; %atmospheric pressure in N/m^2 (value used in Sidebotham)
R_specific = 287.058;  %J/kg/K for dry air
a = 0.25;   %fraction of solar radiation from sun reflected back by earth (albedo) (Cooper1985)

%%%%%%%%%%% Bee Parameters same for BB and HB %%%%%%%%%%%%%%%%
s = 0.9965;  %ratio calculated from Church1960 data converted to K
c = 3.349;  %specific heat (0.8 cal/g*degC converted to J/g*degC *4.1868), cited in May1976
C_l = 2.429809*10^(-7);   %fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
n = 1.975485;       %%fitted from log(Nu) = log(Re), or Nu = C_le^n with CChurch1960 data
delta_T_h = 2.9;
alpha_si = 0.25;     %shape factor for incoming solar radiation (Cooper1985)
alpha_so = 0.5;     %fraction of surface of bee that is irradiated with outgoing solar radiation (Cooper1985)
alpha_th = 0.5;     %fraction of surface of bee that is irradiated with thermal radiation (Cooper1985)
r = 0.0367/9;  %Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
R_0 = 0.000616/2;  %radius of nectar droplet, half avg width of tongue, in m, so drop is width of tongue

D_fl = 0.95;  %percent of distance between resting and flying for metabolic rate and E


if Bumblebee==true
    A_th = 9.3896*10^(-5)  ; %thorax surface area in m^2, from Church1960
    A_h = 3.61375*10^(-5)  ; %head surface area in m^2, own data
    M_b = 0.149;   %mass of the bee in g, Joos1991, default
    M_th = 0.057;  %mass of thorax in g, Joos1991
    l_th = 0.005467;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    epsilon_a = 0.935;   % absorptivity of bees (Willmer1981, ,te)
    v_options = [0 0.1 4.1 2.05 1];   %make the bee be out of wind when resting/shivering, default
    epsilon_e = 0.97;       %(fill in the reference for this!)
    T_mK = 42+273.15;     %median temp for  abdomen cooling
    I_resting = 0.001349728;     %Kammer1974, table 1, for 25C, converted to W
    I_flying = 0.06404973;   %fitted
    I_shivering = I_flying;
%     I_flowers = mean([I_flying,I_resting]);
%     I_flowers = 0.0578;  %90% of distance between values
%     I_flowers = 0.0076;  %10% of distance between values
    I_flowers = D_fl*(I_flying-I_resting)+I_resting;  %fraction of distance between resting and flying that flowers is at
    masses = [0.177 0.177 0.177 0.177 0.177];   %reference weight for Kammer (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15, 25+273.15, 25+273.15];   %Reference temp is 25C for Kammer (resting & flying)
    norm_constants = [I_resting, I_shivering, I_flying, I_flying, I_flowers];   %resting/shivering/flying = 1,2,3
    E_resting = 0.63; 
    E_flying = 0.009176471;
    E_shivering = E_flying;
%     E_flowers = mean([E_resting, E_flying]);
%     E_flowers = 0.0713;  %corresponds to I_flowers at 90% 
%     E_flowers = 0.5679;  %corresponds to I_flowers at 10% 
    E_flowers = (1-D_fl)*(E_resting-E_flying)+E_flying; %solve for E_flowers
    E_opts = [E_resting E_shivering E_flying E_flying E_flowers]; %fitted values

    y0_default = 30+273.15;   %initial temperature of the bee's thorax in K 
    LethalTemp = 45;  %lethal thorax/air temp in C
    CoolingTemp = 42;   %thorax temp where cooling begins
    FlyingTemp = 30;      %thorax temp where flight can begin
    
    if VaryMass==true  %override to scale body and change fight speed
        M_th_default = 0.057;  %mass of thorax in g, Joos1991
        M_b_default = 0.149;   %mass of the bee in g, Joos1991, default
        scale_m = m_thorax/M_th_default;    %how much scaling body mass by 
        l_th = 0.005467*(scale_m^(1/3));   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
        A_th = 4*pi*(l_th/2)^2  ; %thorax surface area in m^2, from Church1960
        A_h = 3.61375*10^(-5)*scale_m^(2/3)  ; %head surface area in m^2, own data
        M_th = m_thorax;  %mass of thorax in g, Joos1991
        M_b = M_b_default*scale_m;   %mass of the bee in g, Joos1991, default
        flying_vs = [1, 4.5, 5.5];
        v_options = [0 0.1 flying_vs(v_bee)];   %make the bee be out of wind when resting/shivering, default
    end
    
    if oilseed==true
        t_handling = 5.03;   %flower handling time
        t_flybetween = 2.23;   %time to fly between flowers
    end
    if fieldbeans==true
        t_handling = 3.4;   %flower handling time
        t_flybetween = 2.23;   %time to fly between flowers
    end

end

if Honeybee==true
   l_th = 0.004;   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
    A_th = 4.5*10^(-5)  ; %thorax surface area in m^2, from ???
    A_h = 2.46*10^(-5)  ; %head surface area in m^2, from Cooper1985 
    M_b = 0.100;   %mass of the bee in g, Joos1991
    M_th = 0.0407;  %mass of thorax in g, Joos1991
    epsilon_a = 0.91;   % absorptivity of bees (Willmer1981, ,te)
    v_options = [0 0.1 5.6 2.8 0.1];   %make the bee be out of wind when resting/shivering, default
    epsilon_e = 0.97;       %(fill in the reference for this!)
    T_mK = 47.9+273.15;     %median temp for  evaporative cooling

    I_resting = 5.65*(80/1000)*(1/1000); %Rothe1989, mW/g -> W, 80mg reference mass
    I_flying = 0.034027;     %%fitted value
%     I_flowers = mean([I_flying,I_resting]);   %fitted
%     I_flowers = 0.0307;  %90% of flying
%     I_flowers = 0.0038;  %10% of distance between resting and flying
    I_flowers = D_fl*(I_flying-I_resting)+I_resting;  %fraction of distance between resting and flying that flowers is at
    masses = [0.08 0.08 0.08 0.08 0.08];   %reference weight for Rothe/Nachtigal (flying)
    RefTemps = [25+273.15, 25+273.15, 25+273.15 25+273.15, 25+273.15];   %Reference temp 
    norm_constants = [I_resting, I_flying, I_flying, I_flying, I_flowers];   %resting/shivering/flying = 1,2,3
    E_resting = 0.63;
    E_flying = 0.008;
    E_shivering = E_flying;
    E_flybetween = E_flying;
%     E_flowers = mean([E_resting, E_flying]);
%     E_flowers = 0.0702;  %corresponds to I_flowers at 90% 
%     E_flowers = 0.5678;  %corresponds to I_flowers at 10% 
    E_flowers = (1-D_fl)*(E_resting-E_flying)+E_flying; %solve for E_flowers
    E_opts = [E_resting E_shivering E_flying E_flybetween E_flowers]; %fitted values

    y0_default = 39+273.15;   %initial temperature of the bee's thorax in K (fill in ref)
    LethalTemp = 52;  %lethal thorax/air temp in C
    CoolingTemp = 47.9;   %thorax temp where cooling begins
    FlyingTemp = 35;      %thorax temp where flight can begin
    
    if VaryMass==true  %override to scale body and change flight speed
        M_th_default = 0.0407;  %mass of thorax in g, Joos1991
        M_b_default = 0.100;   %mass of the bee in g, Joos1991
        scale_m = m_thorax/M_th_default;    %how much scaling body mass by 
        M_th = m_thorax;
        M_b = M_b_default*scale_m;   %mass of the bee in g, Joos1991, default
        l_th = 0.004*(scale_m^(1/3));   %characteristic dimension of thorax in m (avg thorax diam, from Mitchell1976/Church1960)
        A_th = 4*pi*(l_th/2)^2  ; %thorax surface area in m^2, from  Cooper1985
        A_h = 2.46*10^(-5)*scale_m^(2/3)  ; %head surface area in m^2, from Cooper1985 
        flying_vs = [3.6, 4.5, 9];
        v_options = [0 0.1 flying_vs(v_bee)];   %make the bee be out of wind when resting/shivering, default
    end

    if oilseed==true
        t_handling = 4;   %flower handling time
        t_flybetween = 2.23;   %time to fly between flowers
    end
    if fieldbeans==true
        t_handling = 11.9;   %flower handling time
        t_flybetween = 2.23;   %time to fly between flowers
    end

end

%%%%%% Forage Patch values (needed in Set_Rates)
v_flight = v_options(3);
distance_to_patch = 100;   %distance to patch in m
t_bout_minutes = 60;  %60 min bout default for now
t_bout = t_bout_minutes*60;  %transform to seconds
T_0 = FlyingTemp;  %initial thorax temp is the minimum necessary for flight
t_0 = distance_to_patch/v_flight;   %time to fly to field

t_forage = t_bout - 2*t_0;   %time available for foraging
flowervisits_max = floor((t_forage-t_handling)/(t_handling + t_flybetween));   %number of flowers possible to visit
t_bout_actual = 2*t_0+flowervisits_max*(t_handling+t_flybetween)+t_handling;
%using t_forage-t_handling because flower visiting has to start and end
%with a flower visit