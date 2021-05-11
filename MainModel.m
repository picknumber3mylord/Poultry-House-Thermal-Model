% Main Model for NZE Poultry House
% Nathan Shang, Roxy Wilcox, Fermin Banuelos-Gonzalez
% Edited 5/2/2021

clc
clear all

% Thermal model for calculating energy consumption of a model poultry house
% Uses metric units
% Outputs: Thermal energy needed in J
% inputs needed: 
    % climate data
        % average temperature 
        % solar radiation
    % poultry house geometry    

% input reading
inpData2019 = readtable('Thermal_Model_Input_Data_2019.csv','PreserveVariableNames',true);
inpData2018 = readtable('Thermal_Model_Input_Data_2018.csv','PreserveVariableNames',true);
inpData2017 = readtable('Thermal_Model_Input_Data_2017.csv','PreserveVariableNames',true);
inpData2016 = readtable('Thermal_Model_Input_Data_2016.csv','PreserveVariableNames',true);
designData = xlsread('Mk1.xlsx');
% organize input data into matrix

% tData = inpData(%indexing)  % vector of time values % might be unneeded
tempData2019 = table2array(inpData2019(3:end,21))';   % vectors of hourly temperature data from NSRDB
tempData2018 = table2array(inpData2018(3:end,21))';
tempData2017 = table2array(inpData2017(3:end,21))';
tempData2016 = table2array(inpData2016(3:end,21))';

% solar radiation data from NSRDB
DHI2019 = table2array(inpData2019(3:end,6))';    
DHI2018 = table2array(inpData2018(3:end,6))';
DHI2017 = table2array(inpData2017(3:end,6))';
DHI2016 = table2array(inpData2016(3:end,6))';

DNI2019 = table2array(inpData2019(3:end,7))';
DNI2018 = table2array(inpData2018(3:end,7))';
DNI2017 = table2array(inpData2017(3:end,7))';
DNI2016 = table2array(inpData2016(3:end,7))';

GHI2019 = table2array(inpData2019(3:end,8))';
GHI2018 = table2array(inpData2018(3:end,8))';
GHI2017 = table2array(inpData2017(3:end,8))';
GHI2016 = table2array(inpData2016(3:end,8))';

SZA2019 = table2array(inpData2019(3:end,14))';
SZA2018 = table2array(inpData2018(3:end,14))';
SZA2017 = table2array(inpData2017(3:end,14))';
SZA2016 = table2array(inpData2016(3:end,14))';

% extract day numbers
day = table2array(inpData2019(3:end,3));

% define constants necessary for thermal model
rSideInsu = designData(3,1);  %R-value of wall insulation (m^2*K/W)
sideMetalCond = designData(1,1);  %thermal conductivity of metal siding material (W/m/K)
sideMetalThick = designData(5,3);  %thickness of metal siding (m)
aSide =  designData(5,4);  %area of wall minus area of windows (m^2)

Tset = 21; %inside temperature (C)
Tout = [tempData2016; tempData2017; tempData2018; tempData2019];  %outside temperature (C)

CpAir = 1005;  % specific heat of air in (J/kg/K);
rhoAir = 1.204;  % density of air in (kg/m^3); 

pFloor = designData(5,5); % perimeter of flooring (m)
A_f = ; % Area of flooring (m^2)
fFloor = 1.03; %Perimeter heat loss factor (W/m/K)

long = designData(5,1);  %total poultry house length (m)
wide = designData(5,2);  %total poultry house width (m)

rRoofInsu = designData(3,2); %R-value of roof insulation (m^2*K/W)
roofMetalCond = sideMetalCond; %thermal conductivity of metal roofing material (W/m/K)
roofMetalThick = sideMetalThick; %thickness of metal roofing (m)
aRoof = designData(5,6);%area of roof (m^2)

h_is = 3.45; % heat transfer coefficient between air and surface (W/m^2/K)
A_at = 4.5; % dimensionless ratio b/t internal surfaces area and floor area
            % assumed value in accordance with ISO 17390

sensiDay = 4.1; %Senisible heat output of birds during lit hours (W/kg)
sensiNight = 3.2; %sensible heat output of birds during dark hours (W/kg)
chickWeight = 1.713; %weight of chickens (kg)
numChicken = 69376; %Number of chickens in the house
lightOn = 5; %time in hours that house lights are turned on
lightOff = lightOn + 15; %time in hours that house lights are turned off

ventTilt = 90; %tilt of vents with respect to horizontal in degrees
longitude = 121.74; %longitude of Davis, CA in decimal degrees West
latitude = 38.53; %latitude of Davis, CA in decimal degrees
deltaTutc = 8; %Difference between PST and UTC
GHI = [GHI2016; GHI2017; GHI2018; GHI2019];
SZA = [SZA2016; SZA2017; SZA2018; SZA2019];
SHGCvent = 0.13; %Solar Heat Gain Coefficient of Louver vents
aVent = designData(5,7); %area of vents (m^2)

%Energy Consumers%
eggPackerE = 2.5; %Energy consumed by egg packer (kW*hr/day)
eggBeltE = 237.44; %Energy consumed by egg belt (kW*hr/day)
eggElevatorE = 3.6; %Energy consumed by egg elevator (kW*hr/day)
eggWasherE = 175; %Energy consumed by egg washer (kW*hr/day)
lightsE = 13.399; %Energy consumed by lighting (kW*hr/day)
manureBeltE = 20.768; %Energy consumed by manure belt (kW*hr/day)
feedingBeltE = 54.83 %Energy consumed by feeding belt (kWhr/day)

%below is basic calculations to do that are used in functions but remain
%constant, no matter the time
rSideMetal = sideMetalThick/sideMetalCond;
uSide = 1/(rSideInsu+rSideMetal);
ventRate = 12*numChicken/3600;  % ventilation rate (m^3/s)
machinaryE = (eggPackerE + eggBeltE + eggElevatorE + eggWasherE + lightsE + manureBeltE + feedingBeltE) * 1000; %Sums all machinary energy costs and converts them into W*hrs/day

rRoofMetal = roofMetalThick/roofMetalCond; %Calculates R-value of roof metal (m^2*K/W)
uRoof = 1/(rRoofMetal + rRoofInsu); %Overall heat transfer coefficient of room (W/m^2/K)

tData = [0:length(tempData2019)-1];

energyCost2019 = calcCost(tempData2019);
energyCost2018 = calcCost(tempData2018);
energyCost2017 = calcCost(tempData2017);
energyCost2016 = calcCost(tempData2016);
monthlyUse([energyCost2016; energyCost2017; energyCost2018; energyCost2019]);
yearlyUse([energyCost2016; energyCost2017; energyCost2018; energyCost2019]);

% plotting
% plotting hourly data
figure(1)
plot(tData, energyCost2019, 'r', tData, energyCost2018, 'g', tData, energyCost2017, 'b', tData, energyCost2016, 'k')
legend('2019', '2018', '2017', '2016')
title('Hourly Energy Use from 2016-2019')
xlabel('hours passed')
ylabel('Energy Usage (kWh)')

% function for calculating energy requirements
% employs method described by ISO 13790
% modified for poultry house by Costantino et al.
% heat transfer coefficient represent with "H" (W/K)
% T represents temperature in C
% output: energy used to heat/ventilate the house in W
function energyUse = calcCost(tempData)

    energyCost = [];
    mE = machineryE / 24; % converts machinery energy consumption to Whr/hr
    
    for i = 1:length(tempData)
        heatNeed = ;%however we are going to check temperature
        if heatNeed <= 0
            energyCost(i) = mE + fanE + evapCoolE% + whatever else is needed;
        else
            heatCost = (1/0.92) * heatNeed;
            energyCost(i) = heatCost + mE + fanE
        end
            
end

%function for finding heat loss through walls
%inputs: overall heat transfer coeff of walls, area of wall, temperature inside, temperature outside
%output: energy loss in Watts
function wallEnergy = wallE(uSide, aSide, Tset, Tout)
    wallEnergy = uSide*aSide*(Tset-Tout); %energy loss through walls
end

% function for finding heat loss for ventilation
% inputs: specific heat of air, density of air
          % venitlation rate, temperature inside/outside
% output: energy loss in Watts
function ventEnergy = ventE(Cp, rhoAir, ventRate, Tset, Tout)
    ventEnergy = Cp*rhoAir*ventRate*(Tset - Tout); % energy loss to ventilation
end

%function for finding energy loss through floor
%inputs: perimeter of floor (m), F value of floor (W/m/F), temperature inside (C), temperature outside (C)
%outputs: energy loss in Watts
 function floorEnergy = floorE(pFloor, fFloor, Tset, Tout)
    floorEnergy = pFloor*fFloor*(Tset-Tout); %energy loss through floor
 end

 %function for finding energy loss through roof
 %inputs: overall heat transfer coefficient of room, roof area, temperature
        %inside, temperature outside
 %outputs: energy loss in Watts
 function roofEnergy = roofE(uRoof, aRoof, Tset, Tout)
    roofEnergy = uRoof*aRoof*(Tset - Tout); %energy loss through roof
 end
 
%function for finding heat production from chickens
%inputs: current time, sensible heat production of birds during day, 
        %during night, weight of birds, time lights on, time lights off
%outputs: energy production in Watts
function chickEnergy = chickE(t, sensiDay, sensiNight, chickWeight, lightOn, lightOff, numChicken)
    time = mod(t,24);
    if ((time >= lightOn) && (time < lightOff))
        chickEnergy = numChicken*sensiDay*chickWeight; %during lit hours, use senisble heat production for daytime
    else 
        chickEnergy = numChicken*sensiNight*chickWeight; %during dark hours, use sensible heat production for night
    end 
end   

%function for calculating solar radiation energy gain 
%inputs: SHGC, area of surface, tilt of surface, longitude, latitude, day of year, time, GHI, difference between local time and UTC time, zenith angle of sun
%outputs: Energy gain in Watts
function solarEnergy = solarE(SHGC, area, tilt, longitude, latitude, day, t, GHI, deltaTutc, sunZenithAngle)
    if (sunZenithAngle == 90)
        sunZenithAngle = 91; %% sun zenith angle of 90 means the sun is setting and confuses the code. Make it 91 to avoid breaking it
    end
    solarRadS = solarR(tilt, longitude, latitude, day, t, GHI, deltaTutc, 0, sunZenithAngle); %total solar radiation coming from the South in W/m^2
    solarRadN = solarR(tilt, longitude, latitude, day, t, GHI, deltaTutc, 180, sunZenithAngle); %total solar radiation coming from the North in W/m^2
    if (solarRadS<0)
        solarRadS = 0; %Because we cannot have a negative solar radiation, these if statements act as a check
    end
    if (solarRadN<0)
       solarRadN = 0; 
    end
    solarEnergyS = SHGC*solarRadS*area/2; %energy gain from South solar radiation in W
    solarEnergyN = SHGC*solarRadN*area/2; %energy gain from North solar radiation in W
    solarEnergy = solarEnergyS + solarEnergyN; %energy gain from solar radiation in W
end

%function for calculating solar radiation incident on surface
%inputs: tilt angle in degrees, longitude, latitude, day of year, time Matrix, GHI, difference between local time and UTC time, Azimuth angle of module, zenith angle of sun
%outputs: Solar radiation incident on surface in W/m^2
function solarRad = solarR(tilt, longitude, latitude, day, t, GHI, deltaTutc, surfaceAzimuthAngle, sunZenithAngle)
    time = (mod(t,24) + 0.5); %Get hourly times and add a half hour because NSRDB data collected on the half hour
    B = (day-1)*360/365; %Calculation of constant for future use
    Lst = deltaTutc*15; %Calculation of standard meridian for time zone
    E = 229.2*(0.000075 + 0.001868*cosd(B) - 0.032077*sind(B) - 0.014615*cosd(2*B) - 0.04089*sind(2*B));
    timeChangeMin = 4*(Lst - longitude) + E; %Caclulation of change in time from regular to solar in minutes
    Tsol = time + timeChangeMin/60; %Calculation of solar time
    w = 15*(Tsol-12); %Calculation of hour angle
    delta = (180/pi)*(0.006918 - 0.399912*cosd(B) + 0.070257*sind(B) - 0.006758*cosd(2*B) + 0.000907*sind(2*B) - 0.002697*cosd(3*B) + 0.00148*sind(3*B)); %Calculation of declination angle
    theta = acosd(sind(delta)*sind(latitude)*cosd(tilt) - sind(delta)*cosd(latitude)*sind(tilt)*cosd(surfaceAzimuthAngle) + cosd(delta)*cosd(latitude)*cosd(tilt)*cosd(w) + cosd(delta)*sind(latitude)*sind(tilt)*cosd(surfaceAzimuthAngle)*cosd(w) + cosd(delta)*sind(tilt)*sind(surfaceAzimuthAngle)*sind(w)); %calculation of angle of incidence
    solarRad = GHI*(cosd(theta)/cosd(sunZenithAngle));%calculation of amount of solar radiation incident on surface
end

function monthlyUse(enrgCost)
    enrgCost = enrgCost / 1000; % converted to kWh
    for i = 1:4
        jan = sum(enrgCost(i:end,1:(24*31));
        feb = sum(enrgCost(i:end,(24*31+1):(24*59));
        mar = sum(enrgCost(i:end,(24*59+1):(24*90));
        apr = sum(enrgCost(i:end,(24*90+1):(24*120));
        may = sum(enrgCost(i:end,(24*120+1):(24*151));
        jun = sum(enrgCost(i:end,(24*151+1):(24*181));
        jul = sum(enrgCost(i:end,(24*181+1):(24*212));
        aug = sum(enrgCost(i:end,(24*212+1):(24*243));
        sep = sum(enrgCost(i:end,(24*243+1):(24*273));
        oct = sum(enrgCost(i:end,(24*273+1):(24*304));
        nov = sum(enrgCost(i:end,(24*304+1):(24*334));
        dec = sum(enrgCost(i:end,(24*334+1):(24*365));
        fprintf('The average energy usage in January was %d kWh.\n', jan, 'all')/4)
        fprintf('The average energy usage in February was %d kWh.\n', feb, 'all')/4)
        fprintf('The average energy usage in March was %d kWh.\n', mar, 'all')/4)
        fprintf('The average energy usage in April was %d kWh.\n', apr, 'all')/4)
        fprintf('The average energy usage in May was %d kWh.\n', may, 'all')/4)
        fprintf('The average energy usage in June was %d kWh.\n', jun, 'all')/4)
        fprintf('The average energy usage in July was %d kWh.\n', jul, 'all')/4)
        fprintf('The average energy usage in August was %d kWh.\n', aug, 'all')/4)
        fprintf('The average energy usage in September was %d kWh.\n', sep, 'all')/4)
        fprintf('The average energy usage in October was %d kWh.\n', oct, 'all')/4)
        fprintf('The average energy usage in November was %d kWh.\n', nov, 'all')/4)
        fprintf('The average energy usage in December was %d kWh.\n', dec, 'all')/4)
        figure(i+1)
        bar([jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec])
        title(i)
end
    
function yearlyUse(enrgCost)
    enrgCost = enrgCost / 1000; % convert to kWh
    perYearUse = sum(enrgCost, 2);
    fprintf('The energy usage in 2016 was %d kWh.\n', perYearUse(1));
    fprintf('The energy usage in 2017 was %d kWh.\n', perYearUse(2));
    fprintf('The energy usage in 2018 was %d kWh.\n', perYearUse(3));
    fprintf('The energy usage in 2019 was %d kWh.\n', perYearUse(4));
end
