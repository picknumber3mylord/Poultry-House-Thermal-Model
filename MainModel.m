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
rho = %density of air (kg/m^3)
ventRate = %ventilation rate (m^3/s)
Tset = %inside temperature (C)
Tout = %outside temperature (C)

% function for finding heat loss for ventilation
% inputs: specific heat of air, density of air
          % venitlation rate, temperature inside/outside
% output: energy loss in Joules 
function ventEnergy = ventE(Cp, rho, ventRate, Tset, Tout)
    ventEnergy = Cp*rho*ventRate*(Tset - Tout); % energy consumption
end

% function for converting Joules to kWh
function kWH = convKWH(inp)
    kWH = inp * 3.6 * 10^6;
end





