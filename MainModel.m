% Main Model for NZE Poultry House
% Nathan Shang, Roxy Wilcox, Fermin Banuelos-Gonzalez
% Edited 2/23/2021

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
inpData = xlsx
% organize input data into matrix

% define constants necessary for thermal model

Cp = %specific heat of air (J/kg/K)
ventRate = %ventilation rate (m^3/s)
Tset = %inside temperature (C)
Tout = %outside temperature (C)
Tearth = %temperature of ground (C)

floorThick = %floor thickness (m)
floorCond = %thermal conductivity of flooring matieral (W/m/K)
long = %total poultry house length (m)
wide = %total poultry house width (m)

rRoofInsu = %R-value of roof insulation (m^2*K/W)
roofMetalCond = %thermal conductivity of metal roofing material (W/m/K)
roofMetalThick = %thickness of metal roofing (m)
aRoof = %area of roof (m^2)

sensiDay = %Senisible heat output of birds during lit hours (W/kg)
sensiNight = %sensible heat output of birds during dark hours (W/kg)
chickWeight = %weight of chickens (kg)
litTime = %number of hours house is lit for
lightOn = %time in hours that house lights are turned on
lightOff = %time in hours that house lights are turned off


%below is basic calculations to do that are used in functions but remain
%constant, no matter the time

rFloor = floorThick/floorCond; %Calculates R-value of given floor (m^2*K/W)
uFloor = 1/rFloor; %Overall heat transfer coefficient (W/m^2/K)
aFloor = wide*long; %Calculates area of floor (m^2)

rRoofMetal = roofMetalThick/roofMetalCond; %Calculates R-value of roof metal (m^2*K/W)
uRoof = 1/(rRoofMetal + rRoofInsu); %Overall heat transfer coefficient of room (W/m^2/K)

% function for finding heat loss for ventilation
% inputs: specific heat of air, density of air
          % venitlation rate, temperature inside/outside
% output: energy loss in Joules 
function ventEnergy = ventE(Cp, ventRate, Tset, Tout)
    ventEnergy = Cp*densCalc(Tset)*ventRate*(Tset - Tout); % energy loss to ventilation
end

%function for finding energy loss through floor
%inputs: overall heat transfer coefficient of floor, area of floor, temperature inside, 
        %temperature of ground
%outputs: energy loss in Joules
 function floorEnergy = floorE(uFloor, aFloor, Tset, Tearth)
    floorEnergy = aFloor*uFloor*(Tset-Tearth); %energy loss through floor
 end

 %function for finding energy loss through roof
 %inputs: overall heat transfer coefficient of room, roof area, temperature
        %inside, temperature outside
 %outputs: energy loss in Joules
 function roofEnergy = roofE(uRoof, aRoof, Tset, Tout)
    roofEnergy = uRoof*aRoof*(Tset - Tout); %energy loss through roof
 end
 
%function for finding heat production from chickens
%inputs: current time, sensible heat production of birds during day, 
        %during night, weight of birds, time lights on, time lights off
%outputs: energy production in Joules
function chickEnergy = chickE(time, sensiDay, sensiNight, chickWeight, litTime, lightOn, lightOff)
    litTime = litTime*60*60; %Converts litTime from Hours to Seconds
    if ((time >= lightOn) && (time < lightOff))
        chickEnergy = sensiDay*chickWeight*litTime; %during lit hours, use senisble heat production for daytime, outputs in J
    else 
        chickEnergy = sensiNight*chickWeight*litTime; %during dark hours, use sensible heat production for night, outputs in J
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

% function for converting Joules to kWh
function kWH = convKWH(inp)
    kWH = inp * 3.6 * 10^6;
end





