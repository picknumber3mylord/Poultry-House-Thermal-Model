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

CpAir = 1000 (J/kg/K);
rhoAir = 1.225 (kg/m^3);
ventRate = inpData(%indexing);  % ventilation rate (m^3/s)

% Don't think this will be needed
% main loop will give inputs from tempData 
% Tset = % inside temperature (C)
% Tout = % outside temperature (C)

% looping through NSRDB data points
energyCost2019 = 0;
energyCost2018 = 0;
energyCost2017 = 0;
energyCost2016 = 0;
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





