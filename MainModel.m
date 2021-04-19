% Main Model for NZE Poultry House
% Nathan Shang, Roxy Wilcox, Fermin Banuelos-Gonzalez
% Edited 4/19/2021

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
fFloor = 0.72; %F value of floor from ASHRAE 90.1 (BTU/hr/ft/F)

long = designData(5,1);  %total poultry house length (m)
wide = designData(5,2);  %total poultry house width (m)

rRoofInsu = designData(3,2); %R-value of roof insulation (m^2*K/W)
roofMetalCond = sideMetalCond; %thermal conductivity of metal roofing material (W/m/K)
roofMetalThick = sideMetalThick; %thickness of metal roofing (m)
aRoof = designData(5,6);%area of roof (m^2)

sensiDay = 4.1; %Senisible heat output of birds during lit hours (W/kg)
sensiNight = 3.2; %sensible heat output of birds during dark hours (W/kg)
chickWeight = 1.713; %weight of chickens (kg)
numChicken = 69376; %Number of chickens in the house
lightOn = 5; %time in hours that house lights are turned on
lightOff = lightOn + 15; %time in hours that house lights are turned off

ventTilt = 90; %tilt of vents with respect to horizontal in degrees
latitude = 38.5449; %latitude of Davis, CA in decimal degrees
GHI = [GHI2016; GHI2017; GHI2018; GHI2019];
SHGCvent = 0.13; %Solar Heat Gain Coefficient of Louver vents
aVent = designData(5,7); %area of vents (m^2)

%Energy Consumers%
eggPackerE = 2.5; %Energy consumped by egg packer (kW*hr/day)
eggBeltE = 237.44; %Energy consumped by egg belt (kW*hr/day)
eggElevatorE = 3.6; %Energy consumped by egg elevator (kW*hr/day)
eggWasherE = 175; %Energy consumped by egg washer (kW*hr/day)
lightsE = 13.399; %Energy consumped by lighting (kW*hr/day)
manureBeltE = 20.768; %Energy consumped by manure belt (kW*hr/day)

%below is basic calculations to do that are used in functions but remain
%constant, no matter the time
rSideMetal = sideMetalThick/sideMetalCond;
uSide = 1/(rSideInsu+rSideMetal);
ventRate = 12*numChicken/3600;  % ventilation rate (m^3/s)
machinaryE = (eggPackerE + eggBeltE + eggElevatorE + eggWasherE + lightsE + manureBeltE) * 1000; %Sums all machinary energy costs and converts them into W*hrs/day

fFloor = fFloor*1055.06*3.28084/3600; %Coverts F value of floor from BTU/hr/ft/F to W/m/F
%1055.06 J = 1 BTU, 3.28084 ft = 1 m; 1 hr = 3600 sec

rRoofMetal = roofMetalThick/roofMetalCond; %Calculates R-value of roof metal (m^2*K/W)
uRoof = 1/(rRoofMetal + rRoofInsu); %Overall heat transfer coefficient of room (W/m^2/K)

% looping through NSRDB data points
energyCost2019 = [];
energyCost2018 = [];
energyCost2017 = [];
energyCost2016 = [];
tData = [0:length(tempData2019)-1];

for i = 2:length(tData)
    % do calcs
    energyCost2016(i) = wallE(uSide, aSide, Tset, Tout(1,i)) + ventE(CpAir, rhoAir, ventRate, Tset, Tout(1,i)) + floorE(pFloor, fFloor, Tset, Tout(1,i)) + roofE(uRoof, aRoof, Tset, Tout(1,i)) + chickE(tData(i), sensiDay, sensiNight, chickWeight, lightOn, lightOff, numChicken) + solarE(SHGCvent, ventTilt, latitude, day(i), GHI(1,i), aVent);
                    
    energyCost2017(i) = wallE(uSide, aSide, Tset, Tout(2,i)) + ventE(CpAir, rhoAir, ventRate, Tset, Tout(2,i)) + floorE(pFloor, fFloor, Tset, Tout(2,i)) + roofE(uRoof, aRoof, Tset, Tout(2,i)) + chickE(tData(i), sensiDay, sensiNight, chickWeight, lightOn, lightOff, numChicken) + solarE(SHGCvent, ventTilt, latitude, day(i), GHI(2,i), aVent);
                    
    energyCost2018(i) = wallE(uSide, aSide, Tset, Tout(3,i)) + ventE(CpAir, rhoAir, ventRate, Tset, Tout(3,i)) + floorE(pFloor, fFloor, Tset, Tout(3,i)) + roofE(uRoof, aRoof, Tset, Tout(3,i)) + chickE(tData(i), sensiDay, sensiNight, chickWeight, lightOn, lightOff, numChicken) + solarE(SHGCvent, ventTilt, latitude, day(i), GHI(3,i), aVent);
                    
    energyCost2019(i) = wallE(uSide, aSide, Tset, Tout(4,i)) + ventE(CpAir, rhoAir, ventRate, Tset, Tout(4,i)) + floorE(pFloor, fFloor, Tset, Tout(4,i)) + roofE(uRoof, aRoof, Tset, Tout(4,i)) + chickE(tData(i), sensiDay, sensiNight, chickWeight, lightOn, lightOff, numChicken) + solarE(SHGCvent, ventTilt, latitude, day(i), GHI(4,i), aVent);
                    
end

monthlyUse([energyCost2016; energyCost2017; energyCost2018; energyCost2019]);
yearlyUse([energyCost2016; energyCost2017; energyCost2018; energyCost2019]);

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
    TsetF = tempToFahr(Tset); %converts Tset to deg Farhenheit
    ToutF = tempToFahr(Tout); %converts Tout to deg Farhenheit
    tempDiff = TsetF-ToutF; %calculates temperature difference
    floorEnergy = pFloor*fFloor*(tempDiff); %energy loss through floor
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
%inputs: SHGC, tilt of surface, latitude, day of year, GHI, area of surface
%outputs: Energy gain in Watts
function solarEnergy = solarE(SHGC, tilt, latitude, day, GHI, area)
    solarR = solarRad(tilt, latitude, day, GHI); %total solar radiation in W/m^2
    solarEnergy = SHGC*solarR*area; %energy gain from solar radiation in W
end

%function for calculating solar radiation incident on surface
%inputs: tilt angle in degrees, latitude, day of year, GHI
%outputs: Solar radiation incident on surface in W/m^2
function solarR = solarRad(tilt, latitude, day, GHI)
    delta = 23.45*sind((360/365)*(284+day));
    alpha = 90 - latitude + delta;
    solarR = (GHI*sind(alpha+tilt))/sind(alpha);
end

%function to convert from deg Celcius to deg Fahrenheit
%inputs: temperature (C)
%outputs: temperature (F)
function tempF = tempToFahr(tempC)
    tempF = tempC*9/5+32;
end

function monthlyUse(enrgCost)
    enrgCost = enrgCost / 1000; % converted to kWh
    fprintf('The average energy usage in January was %d kWh.\n', sum(enrgCost(1:end,1:(24*31)), 'all')/4)
    fprintf('The average energy usage in February was %d kWh.\n', sum(enrgCost(1:end,(24*31+1):(24*59)), 'all')/4)
    fprintf('The average energy usage in March was %d kWh.\n', sum(enrgCost(1:end,(24*59+1):(24*90)), 'all')/4)
    fprintf('The average energy usage in April was %d kWh.\n', sum(enrgCost(1:end,(24*90+1):(24*120)), 'all')/4)
    fprintf('The average energy usage in May was %d kWh.\n', sum(enrgCost(1:end,(24*120+1):(24*151)), 'all')/4)
    fprintf('The average energy usage in June was %d kWh.\n', sum(enrgCost(1:end,(24*151+1):(24*181)), 'all')/4)
    fprintf('The average energy usage in July was %d kWh.\n', sum(enrgCost(1:end,(24*181+1):(24*212)), 'all')/4)
    fprintf('The average energy usage in August was %d kWh.\n', sum(enrgCost(1:end,(24*212+1):(24*243)), 'all')/4)
    fprintf('The average energy usage in September was %d kWh.\n', sum(enrgCost(1:end,(24*243+1):(24*273)), 'all')/4)
    fprintf('The average energy usage in October was %d kWh.\n', sum(enrgCost(1:end,(24*273+1):(24*304)), 'all')/4)
    fprintf('The average energy usage in November was %d kWh.\n', sum(enrgCost(1:end,(24*304+1):(24*334)), 'all')/4)
    fprintf('The average energy usage in December was %d kWh.\n', sum(enrgCost(1:end,(24*334+1):(24*365)), 'all')/4)
end
    
function yearlyUse(enrgCost)
    enrgCost = enrgCost / 1000; % convert to kWh
    perYearUse = sum(enrgCost, 2);
    fprintf('The energy usage in 2016 was %d kWh.\n', perYearUse(1));
    fprintf('The energy usage in 2017 was %d kWh.\n', perYearUse(2));
    fprintf('The energy usage in 2018 was %d kWh.\n', perYearUse(3));
    fprintf('The energy usage in 2019 was %d kWh.\n', perYearUse(4));
end
