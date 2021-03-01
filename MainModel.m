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
inpData = readtable('Thermal_Model_Input_Data.xlsx')
% organize input data into matrix

tData = inpData(%indexing)  % vector of time values
tempData = inpData(%indexing)   % matrix of hourly temperature data from NSRDB
                                % each row corresponds to year
radData = inpData(%indexing)    % solar radiation data from NSRDB

% define constants necessary for thermal model

CpAir = 1000 (J/kg/K);
rhoAir = 1.225 (kg/m^3);
ventRate = inpData(%indexing);  % ventilation rate (m^3/s)

% Don't think this will be needed
% main loop will give inputs from tempData 
% Tset = % inside temperature (C)
% Tout = % outside temperature (C)

% looping through NSRDB data points
energyCost = 0;
for i = 1:length(tData)
    % do calcs
    energyCost += ventE(CpAir, rhoAir, ventRate, Tset, tempData(%indexing))
    energyCost += everything else
end



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





