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

% function for finding heat loss for ventilation
% inputs: specific heat of air, density of air
          % venitlation rate, temperature inside/outside
% output: energy consumption in W

function energyNeeded = ventHeat(% inputs)
    % define constants in equation
    % do math
    energyNeeded = % energy consumption
end

% function for converting Joules to kWh
function kWH = convKWH(inp)
    kWH = inp * 3.6 * 10^6;
end





