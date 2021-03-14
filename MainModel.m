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
inpData2019 = csvread('Thermal_Model_Input_Data_2019.csv', 4, 5);
inpData2018 = csvread('Thermal_Model_Input_Data_2018.csv', 4, 5);
inpData2017 = csvread('Thermal_Model_Input_Data_2017.csv', 4, 5);
inpData2016 = csvread('Thermal_Model_Input_Data_2016.csv', 4, 5);
designData = xlsread(
% organize input data into matrix

% tData = inpData(%indexing)  % vector of time values % might be unneeded
tempData2019 = inpData2019(:,16)';   % vectors of hourly temperature data from NSRDB
tempData2018 = inpData2018(:,16)';
tempData2017 = inpData2017(:,16)';
tempData2016 = inpData2016(:,16)';

% solar radiation data from NSRDB
DHI2019 = inpData2019(:,1)';    
DHI2018 = inpData2018(:,1)';
DHI2017 = inpData2017(:,1)';
DHI2016 = inpData2016(:,1)';

DNI2019 = inpData2019(:,2)';
DNI2018 = inpData2018(:,2)';
DNI2017 = inpData2017(:,2)';
DNI2016 = inpData2016(:,2)';

GHI2019 = inpData2019(:,3)';
GHI2018 = inpData2018(:,3)';
GHI2017 = inpData2017(:,3)';
GHI2016 = inpData2016(:,3)';

SZA2019 = inpData2019(:,9)';
SZA2018 = inpData2018(:,9)';
SZA2017 = inpData2017(:,9)';
SZA2016 = inpData2016(:,9)';

% define constants necessary for thermal model


%rSideInsu = %R-value of wall insulation (m^2*K/W)
%sideMetalCond = %thermal conductivity of metal siding material (W/m/K)
%sideMetalThick = %thickness of metal siding (m)
%aSide = %area of wall minus area of windows (m^2)

%Cp = %specific heat of air (J/kg/K)
%ventRate = %ventilation rate (m^3/s)
%Tset = %inside temperature (C)
%Tout = %outside temperature (C)
=======
%CpAir = 1000 (J/kg/K);
%rhoAir = 1.225 (kg/m^3);
%ventRate = inpData(%indexing);  % ventilation rate (m^3/s)

% Don't think this will be needed
% main loop will give inputs from tempData 
% Tset = % inside temperature (C)
% Tout = % outside temperature (C)

% looping through NSRDB data points
%energyCost2019 = 0;
%energyCost2018 = 0;
%energyCost2017 = 0;
%energyCost2016 = 0;
%for i = 1:length(tData)
    % do calcs
    %energyCost += ventE(CpAir, rhoAir, ventRate, Tset, tempData(%indexing))
    %energyCost += everything else
%end




pFloor = %perimeter of flooring (m)
fFloor = %F value of floor from ASHRAE 90.1 (BTU/hr/ft/F)

long = %total poultry house length (m)
wide = %total poultry house width (m)

rRoofInsu = %R-value of roof insulation (m^2*K/W)
roofMetalCond = %thermal conductivity of metal roofing material (W/m/K)
roofMetalThick = %thickness of metal roofing (m)
aRoof = %area of roof (m^2)

sensiDay = %Senisible heat output of birds during lit hours (W/kg)
sensiNight = %sensible heat output of birds during dark hours (W/kg)
chickWeight = %weight of chickens (kg)
lightOn = %time in hours that house lights are turned on
lightOff = %time in hours that house lights are turned off


%below is basic calculations to do that are used in functions but remain
%constant, no matter the time
rSideMetal = sideMetalThick/sideMetalCond;
uSide = 1/(rSideInsu+rSideMetal);

fFloor = fFloor*1055.06*3.28084/3600; %Coverts F value of floor from BTU/hr/ft/F to W/m/F
%1055.06 J = 1 BTU, 3.28084 ft = 1 m; 1 hr = 3600 sec

rRoofMetal = roofMetalThick/roofMetalCond; %Calculates R-value of roof metal (m^2*K/W)
uRoof = 1/(rRoofMetal + rRoofInsu); %Overall heat transfer coefficient of room (W/m^2/K)

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
 function floorEnergy = floorE(pFloor, fFloor, Tset, Tearth)
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
function chickEnergy = chickE(time, sensiDay, sensiNight, chickWeight, lightOn, lightOff)
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


% function for converting Joules to kWh
function kWH = convKWH(inp)
    kWH = inp * 3.6 * 10^6;
end





