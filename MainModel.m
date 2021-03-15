% Main Model for NZE Poultry House
% Nathan Shang, Roxy Wilcox, Fermin Banuelos-Gonzalez
% Edited 3/1/2021

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
tempData2019 = table2array(inpData2019(4:end,21))';   % vectors of hourly temperature data from NSRDB
tempData2018 = table2array(inpData2018(4:end,21))';
tempData2017 = table2array(inpData2017(4:end,21))';
tempData2016 = table2array(inpData2016(4:end,21))';

% solar radiation data from NSRDB
DHI2019 = table2array(inpData2019(4:end,6))';    
DHI2018 = table2array(inpData2018(4:end,6))';
DHI2017 = table2array(inpData2017(4:end,6))';
DHI2016 = table2array(inpData2016(4:end,6))';

DNI2019 = table2array(inpData2019(4:end,7))';
DNI2018 = table2array(inpData2018(4:end,7))';
DNI2017 = table2array(inpData2017(4:end,7))';
DNI2016 = table2array(inpData2016(4:end,7))';

GHI2019 = table2array(inpData2019(4:end,8))';
GHI2018 = table2array(inpData2018(4:end,8))';
GHI2017 = table2array(inpData2017(4:end,8))';
GHI2016 = table2array(inpData2016(4:end,8))';

SZA2019 = table2array(inpData2019(4:end,14))';
SZA2018 = table2array(inpData2018(4:end,14))';
SZA2017 = table2array(inpData2017(4:end,14))';
SZA2016 = table2array(inpData2016(4:end,14))';

% define constants necessary for thermal model


rSideInsu = designData(3,1);  %R-value of wall insulation (m^2*K/W)
sideMetalCond = designData(1,1);  %thermal conductivity of metal siding material (W/m/K)
sideMetalThick = designData(5,3);  %thickness of metal siding (m)
aSide =  designData(5,4);  %area of wall minus area of windows (m^2)

Cp = 1005; %specific heat of air (J/kg/K)
ventRate = designData(8,1);  %ventilation rate (m^3/s)
Tset = 70; %inside temperature (C)
Tout = [tempData2016; tempData2017; tempData2018; tempData2019];  %outside temperature (C)

CpAir = 1000;  % specific heat of air in (J/kg/K);
rhoAir = 1.225;  % density of air in (kg/m^3);
ventRate = designData(8,1);  % ventilation rate (m^3/s)

% Don't think this will be needed
% main loop will give inputs from tempData 

pFloor = designData(5,5); % perimeter of flooring (m)
fFloor = 0.2; %F value of floor from ASHRAE 90.1 (BTU/hr/ft/F)

long = designData(5,1);  %total poultry house length (m)
wide = designData(5,2);  %total poultry house width (m)

rRoofInsu = designData(3,2); %R-value of roof insulation (m^2*K/W)
roofMetalCond = sideMetalCond; %thermal conductivity of metal roofing material (W/m/K)
roofMetalThick = sideMetalThick; %thickness of metal roofing (m)
aRoof = designData(5,6);%area of roof (m^2)

sensiDay = 4.1; %Senisible heat output of birds during lit hours (W/kg)
sensiNight = 3.2; %sensible heat output of birds during dark hours (W/kg)
chickWeight = 1.713; %weight of chickens (kg)
lightOn = 4; %time in hours that house lights are turned on
lightOff = lightOn + 16; %time in hours that house lights are turned off


%below is basic calculations to do that are used in functions but remain
%constant, no matter the time
rSideMetal = sideMetalThick/sideMetalCond;
uSide = 1/(rSideInsu+rSideMetal);

% fFloor = fFloor*1055.06*3.28084/3600; %Coverts F value of floor from BTU/hr/ft/F to W/m/F
%1055.06 J = 1 BTU, 3.28084 ft = 1 m; 1 hr = 3600 sec

rRoofMetal = roofMetalThick/roofMetalCond; %Calculates R-value of roof metal (m^2*K/W)
uRoof = 1/(rRoofMetal + rRoofInsu); %Overall heat transfer coefficient of room (W/m^2/K)

% looping through NSRDB data points
energyCost2019 = 0;
energyCost2018 = 0;
energyCost2017 = 0;
energyCost2016 = 0;
tData = [0:length(tempData2019)-1] / 24;
for i = 2:length(tData)
    % do calcs
    energyCost2016 = energyCost2016 + wallE(uSide, aSide, Tset, Tout(1,i)) 
                    + ventE(CpAir, ventRate, Tset, Tout(1,i))
                    + floorE(pFloor, fFloor, Tset, Tout(1,i))
                    + roofE(uRoof, aRoof, Tset, Tout(1,i))
                    + chickE(tData(i), sensiDay, sensiNight, chickWeight, lightOn, lightOff);
                    
    energyCost2017 = energyCost2017 + wallE(uSide, aSide, Tset, Tout(2,i)) 
                    + ventE(CpAir, ventRate, Tset, Tout(2,i))
                    + floorE(pFloor, fFloor, Tset, Tout(2,i))
                    + roofE(uRoof, aRoof, Tset, Tout(2,i))
                    + chickE(tData(i), sensiDay, sensiNight, chickWeight, lightOn, lightOff);
                    
    energyCost2018 = energyCost2018 + wallE(uSide, aSide, Tset, Tout(3,i)) 
                    + ventE(CpAir, ventRate, Tset, Tout(3,i))
                    + floorE(pFloor, fFloor, Tset, Tout(3,i))
                    + roofE(uRoof, aRoof, Tset, Tout(3,i))
                    + chickE(tData(i), sensiDay, sensiNight, chickWeight, lightOn, lightOff);
                    
    energyCost2019 = energyCost2019 + wallE(uSide, aSide, Tset, Tout(4,i)) 
                    + ventE(CpAir, ventRate, Tset, Tout(4,i))
                    + floorE(pFloor, fFloor, Tset, Tout(4,i))
                    + roofE(uRoof, aRoof, Tset, Tout(4,i))
                    + chickE(tData(i), sensiDay, sensiNight, chickWeight, lightOn, lightOff);
                    
end

fprintf('The energy usage in 2016 was %d kWh.\n', convKWH(energyCost2016, len(tData)-1);
fprintf('The energy usage in 2017 was %d kWh.\n', convKWH(energyCost2017, len(tData)-1));
fprintf('The energy usage in 2018 was %d kWh.\n', convKWH(energyCost2018, len(tData)-1));
fprintf('The energy usage in 2019 was %d kWh.\n', convKWH(energyCost2019, len(tData)-1));

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
function ventEnergy = ventE(Cp, ventRate, Tset, Tout)
    ventEnergy = Cp*densCalc(Tset)*ventRate*(Tset - Tout); % energy loss to ventilation
end

%function for finding energy loss through floor
%inputs: perimeter of floor (m), F value of floor (W/m/F), temperature inside (C), temperature outside (C)
%outputs: energy loss in Watts
 function floorEnergy = floorE(pFloor, fFloor, Tset, Tout)
    tempDiff = Tset-Tout; %caluclates difference in temperature
    tempDiff = tempToFahr(tempDiff); %converts difference to deg Farhenheit
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
function chickEnergy = chickE(t, sensiDay, sensiNight, chickWeight, lightOn, lightOff)
    time = mod(t,24);
    if ((time >= lightOn) && (time < lightOff))
        chickEnergy = sensiDay*chickWeight; %during lit hours, use senisble heat production for daytime
    else 
        chickEnergy = sensiNight*chickWeight; %during dark hours, use sensible heat production for night
    end 
end

%function to calculate current density of air given T in deg C
%inputs: Indoor Temperature (C)
%outputs: density of air
function dens = densCalc(T)
    R = 287.05; %Universal Gas Constant (J/kg/K)
    T = T + 273.15; %Covert Temp to Kelvin
    P = 101325; %Absolute Pressure (Pa)
    dens = P/(R*T); %Performs Calculation of Density of Air (kg/m^3)
end

%function to convert from deg Celcius to deg Fahrenheit
%inputs: temperature (C)
%outputs: temperature (F)
function tempF = tempToFahr(tempC)
    tempF = tempC*9/5+32;
end


% function for converting W to kWh
function kWH = convKWH(inp, t)
    kWH = (inp / 1000) * t;
end
