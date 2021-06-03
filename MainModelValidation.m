% Main Model for NZE Poultry House
% Nathan Shang, Roxy Wilcox, Fermin Banuelos-Gonzalez
% Edited 5/22/2021

clc
clear all
close all
% Thermal model for calculating energy consumption of a model poultry house
% Uses metric units
% Outputs: Thermal energy needed in J
% inputs needed: 
    % climate data
        % average temperature 
        % solar radiation
    % poultry house geometry    

% input reading
inpData2010 = readtable('Thermal_Model_Input_Data_2010.csv','PreserveVariableNames',true);
inpData2011 = readtable('Thermal_Model_Input_Data_2011.csv','PreserveVariableNames',true);

designData = xlsread('Mk1Indiana.xlsx');
% organize input data into matrix

% tData = inpData(%indexing)  % vector of time values % might be unneeded
tempData2010 = table2array(inpData2010(3:end,21))';
tempData2011 = table2array(inpData2011(3:end,21))';
Tout = [tempData2010; tempData2011];  %outside temperature (C)

%Relative Humiluty data from NSRDB
RHdata2010 = table2array(inpData2010(3:end,20))';
RHdata2011 = table2array(inpData2011(3:end,20))';
RH = [RHdata2010; RHdata2011];

% solar radiation data from NSRDB
GHI2010 = table2array(inpData2010(3:end,8))';
GHI2011 = table2array(inpData2011(3:end,8))';
SZA2010 = table2array(inpData2010(3:end,14))';
SZA2011 = table2array(inpData2011(3:end,14))';
GHI = [GHI2010;GHI2011];
SZA = [SZA2010;SZA2011];

% extract day numbers
day = table2array(inpData2010(3:end,3));
month = table2array(inpData2010(3:end,2));

%HVAC Constants
heatEff = 0.92; %Efficiency of heaters
evapPadEff = .9; %Efficiency of evaporative cooling pads (%)

% define constants necessary for thermal model
rSideInsu = designData(3,1);  %R-value of wall insulation (m^2*K/W)
sideMetalCond = designData(1,1);  %thermal conductivity of metal siding material (W/m/K)
sideMetalThick = designData(5,3);  %thickness of metal siding (m)
aSide =  designData(5,4);  %area of wall minus area of windows (m^2)

% TsetJune = 28.5;
% TsetJuly = 30;
% TsetAug = 29.5;
% TsetSept = 26;
% TsetNov = 23;
% TsetDec = 22;
% TsetJanFeb = 23;
% TsetMar = 25;
% TsetMay = 26;

TsetJanFeb = 23;
TsetMarMay = 25;
TsetJuneAug = 27;
TsetSeptNov = 23+1;
TsetDec = 22;

CpAir = 1005;  % specific heat of air in (J/kg/K);
rhoAir = 1.204;  % density of air in (kg/m^3); 

pFloor = designData(5,5); % perimeter of flooring (m)
fFloor = 1.03; %Perimeter heat loss factor (W/m/K)

long = designData(5,1);  %total poultry house length (m)
wide = designData(5,2);  %total poultry house width (m)

rRoofInsu = designData(3,2); %R-value of roof insulation (m^2*K/W)
roofMetalCond = sideMetalCond; %thermal conductivity of metal roofing material (W/m/K)
roofMetalThick = sideMetalThick; %thickness of metal roofing (m)
aRoof = designData(5,6);%area of roof (m^2)

sensiDay = 4.1; %Senisible heat output of birds during lit hours (W/kg)
sensiNight = 3.2; %sensible heat output of birds during dark hours (W/kg)
chickWeight = 1.713; %weight of chickens (kg)
numChicken = 50000; %Number of chickens in the house
lightOn = 5.5; %time in hours that house lights are turned on
lightOff = 21.75; %time in hours that house lights are turned off

ventTilt = 90; %tilt of vents with respect to horizontal in degrees
longitude = 86.82; %longitude of location in decimal degrees West
latitude = 39.53; %latitude of location in decimal degrees
deltaTutc = 6; %Difference between CST and UTC
SHGCvent = 0.13; %Solar Heat Gain Coefficient of Louver vents
aVent = designData(5,7); %area of vents (m^2)
ventEff = 25.5; %ventilation efficiency (m^3/hr/W)

%Energy Consumers%
eggPackerE = 2.5; %Energy consumed by egg packer (kW*hr/day)
eggBeltE = 237.44; %Energy consumed by egg belt (kW*hr/day)
eggElevatorE = 3.6; %Energy consumed by egg elevator (kW*hr/day)
lightsE = 13.399; %Energy consumed by lighting (kW*hr/day)
manureBeltE = 20.768; %Energy consumed by manure belt (kW*hr/day)
feedingBeltE = 54.83; %Energy consumed by feeding belt (kWhr/day)

%below is basic calculations to do that are used in functions but remain
%constant, no matter the time
rSideMetal = sideMetalThick/sideMetalCond;
uSide = 1/(rSideInsu+rSideMetal);
ventRatePerHen = 0.6;  % ventilation rate (m^3/hen/hr)
ventRate = ventRatePerHen*numChicken/3600; %Ventilation rate (m^3/s)

manureBeltESummer = 7.78*0.1167; %Energy consumed by manure belt (kW*hr/day)
manureBeltEWinter = 7.78*.25;

rRoofMetal = roofMetalThick/roofMetalCond; %Calculates R-value of roof metal (m^2*K/W)
uRoof = 1/(rRoofMetal + rRoofInsu); %Overall heat transfer coefficient of room (W/m^2/K)

% calculating energy usage
tData = [1:length(tempData2010)];
fans = zeros(2,length(tempData2010));
heat = zeros(2,length(tempData2010));
heatNeed = zeros(2,length(tempData2010));
energyCost = zeros(2,length(tempData2010));

for i = 1:2
    for j = 1:length(tempData2010)

%         if (month(j) ==6) %June
%             Tset = TsetJune;
%             manureBeltE = manureBeltESummer;
%         elseif (month(j) == 7) %july
%             Tset = TsetJuly;
%             manureBeltE = manureBeltESummer;            
%         elseif (month(j) == 8)
%             Tset = TsetAug;
%             manureBeltE = manureBeltESummer;
%         elseif (month(j) == 9)
%             Tset = TsetSept;
%             manureBeltE = manureBeltESummer;
%         elseif (month(j) >= 10) && (month(j) <= 11) %btwn oct and nov
%             Tset = TsetNov;
%             manureBeltE = manureBeltEWinter;
%         elseif (month(j) == 12)
%             manureBeltE = manureBeltEWinter;
%             Tset = TsetDec;
%         elseif (month(j) >=1) && (month(j) <= 2) %btwn jan and feb
%             Tset = TsetJanFeb;
%             manureBeltE = manureBeltEWinter;
%         elseif (month(j) == 3)
%             Tset = TsetMar;
%             manureBeltE = manureBeltESummer;
%         else %may
%             Tset = TsetMay;
%             manureBeltE = manureBeltESummer;
%         end
        if (month(j) >= 1) && (month(j) <= 2)
            Tset = TsetJanFeb;
            manureBeltE = manureBeltEWinter;
        elseif (month(j) >= 3) && (month(j) <= 5)
            Tset = TsetMarMay;
            manureBeltE = manureBeltESummer;
        elseif (month(j) >= 6) && (month(j) <= 8)
            Tset = TsetJuneAug;
            manureBeltE = manureBeltESummer;
        elseif (month(j) >= 9) && (month(j) <= 10)
            Tset = TsetSeptNov;
            manureBeltE = manureBeltEWinter;
        elseif (month(j) == 11)
            Tset = 24;
            manureBeltE = manureBeltEWinter;
        else
            Tset = TsetDec;
            manureBeltE = manureBeltEWinter;
        end
        
        machineryE = (eggPackerE + eggBeltE + eggElevatorE + lightsE + manureBeltE + feedingBeltE); %Sums all machinemry energy costs and converts them into kW*hrs/day
        mE = machineryE / 24; % converts machinery energy consumption to kWhr/hr
        

        heatNeed(i,j) = -wallE(uSide, aSide, Tset, Tout(i,j)) - ventE(CpAir, rhoAir, ventRate, Tset, Tout(i,j)) - floorE(pFloor, fFloor, Tset, Tout(i,j)) - roofE(uRoof, aRoof, Tset, Tout(i,j)) + chickE(tData(j), sensiDay, sensiNight, chickWeight, lightOn, lightOff, numChicken) + solarE(SHGCvent, aVent, ventTilt, longitude, latitude, day(j), tData(j), GHI(i,j), deltaTutc, SZA(i,j));
        if heatNeed(i,j) < 0 % net heat lost to external environment
        % need to calculate heating
        % base ventilation rate used
            energyCost(i,j) = mE + (ventRate/ventEff*3600 - heatNeed(i,j)/heatEff) / 1000;
            fans(i,j) = (ventRate/ventEff*3600) / 1000;
            heat(i,j) = -(heatNeed(i,j)/heatEff) / 1000;
        else % net heat gain from external environment
        % cooling is needed from ventilation and evap cooling
            %if (Tset < (Tout(i,j)+2))
                energyCost(i,j) = mE + evapCool(heatNeed(i,j), RH(i,j), Tout(i,j), evapPadEff, CpAir, rhoAir, ventRate, Tset, ventEff);
                fans(i,j) = evapCool(heatNeed(i,j), RH(i,j), Tout(i,j), evapPadEff, CpAir, rhoAir, ventRate, Tset, ventEff);
%             else
%                 fans(i,j) = Cool(heatNeed(i,j), Tout(i,j), CpAir, rhoAir, ventRate, Tset, ventEff);
%                 energyCost(i,j) = mE + fans(i,j);
%                 if (fans(i,j) > 1000000)
%                     energyCost(i,j) = mE + evapCool(heatNeed(i,j), RH(i,j), Tout(i,j), evapPadEff, CpAir, rhoAir, ventRate, Tset, ventEff);
%                     fans(i,j) = evapCool(heatNeed(i,j), RH(i,j), Tout(i,j), evapPadEff, CpAir, rhoAir, ventRate, Tset, ventEff);
%                 end
%             end
            heat(i,j) = 0;
        end
    end
end

% ISU graphing 
janF = sum(fans(2:end,1:(24*31))/31, 'all');
febF = sum(fans(2:end,(24*31+1):(24*59))/28, 'all');
marF = sum(fans(2:end,(24*59+1):(24*90))/31, 'all');
aprF = sum(fans(2:end,(24*90+1):(24*120))/30, 'all');
mayF = sum(fans(2:end,(24*120+1):(24*151))/31, 'all');
junF = sum(fans(2:end,(24*151+1):(24*181))/30, 'all');
julF = sum(fans(2:end,(24*181+1):(24*212))/31, 'all');
augF = sum(fans(2:end,(24*212+1):(24*243))/31, 'all');
sepF = sum(fans(1:end,(24*243+1):(24*273))/30, 'all');
octF = sum(fans(1:end,(24*273+1):(24*304))/31, 'all');
novF = sum(fans(1:end,(24*304+1):(24*334))/30, 'all');
decF = sum(fans(1:end,(24*334+1):(24*365))/31, 'all');
sepF1 = sum(fans(2:end,(24*243+1):(24*273))/30, 'all');
octF1 = sum(fans(2:end,(24*273+1):(24*304))/31, 'all');
novF1 = sum(fans(2:end,(24*304+1):(24*334))/30, 'all');
% monthly heating usage
janH = sum(heat(2:end,1:(24*31))/31, 'all');
febH = sum(heat(2:end,(24*31+1):(24*59))/28, 'all');
marH = sum(heat(2:end,(24*59+1):(24*90))/31, 'all');
aprH = sum(heat(2:end,(24*90+1):(24*120))/30, 'all');
mayH = sum(heat(2:end,(24*120+1):(24*151))/31, 'all');
junH = sum(heat(2:end,(24*151+1):(24*181))/30, 'all');
julH = sum(heat(2:end,(24*181+1):(24*212))/31, 'all');
augH = sum(heat(2:end,(24*212+1):(24*243))/31, 'all');
sepH = sum(heat(1:end,(24*243+1):(24*273))/30, 'all');
octH = sum(heat(1:end,(24*273+1):(24*304))/31, 'all');
novH = sum(heat(1:end,(24*304+1):(24*334))/30, 'all');
sepH1 = sum(heat(2:end,(24*243+1):(24*273))/30, 'all');
octH1 = sum(heat(2:end,(24*273+1):(24*304))/31, 'all');
novH1 = sum(heat(2:end,(24*304+1):(24*334))/30, 'all');
decH = sum(heat(1:end,(24*334+1):(24*365))/31, 'all');
    
y = [sepF manureBeltE lightsE feedingBeltE;
     octF manureBeltE lightsE feedingBeltE;
     novF manureBeltE lightsE feedingBeltE;
     decF manureBeltE lightsE feedingBeltE;
     janF manureBeltE lightsE feedingBeltE;
     febF manureBeltE lightsE feedingBeltE
     marF manureBeltE lightsE feedingBeltE;
     aprF manureBeltE lightsE feedingBeltE;
     mayF manureBeltE lightsE feedingBeltE;
     junF manureBeltE lightsE feedingBeltE
     julF manureBeltE lightsE feedingBeltE;
     augF manureBeltE lightsE feedingBeltE;
     sepF1 manureBeltE lightsE feedingBeltE;
     octF1 manureBeltE lightsE feedingBeltE;
     novF1 manureBeltE lightsE feedingBeltE];
         
figure(1)
bar(y, 'stacked')
% if i == 1
%     title('Energy Usage for 2010')
% elseif i == 2
%     title('Energy Usage for 2011')
% end
xlabel('Month')
ylabel('Energy Usage (kWh/day)')
set(gca,'xticklabel',{'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov'})
legend('Fans', 'Manure Belt', 'Lighting', 'Feeding Belt', 'Location', 'north')
     
     figure(2)
     bar([sepH octH novH decH janH febH marH aprH mayH junH julH augH novH])
    
    title('Heating Energy')
     
     xlabel('Month')
     ylabel('Energy Usage (kWh)')
     set(gca,'xticklabel',{'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Nov'})



monthlyUse(energyCost);


%function for determining required ventilation when not evaporative cooling
function CoolE = Cool(qSup, Tout, CpAir, rhoAir, ventRate, Tset, ventEff)
    qSupNoVent = qSup + ventE(CpAir, rhoAir, ventRate, Tset, Tout); %removal of ventilation from supplemental heating required because it is about to be recalculated (W)
    VentRateReq = qSupNoVent/(CpAir*rhoAir*(Tset - Tout)); %Ventilation rate required given supplemental heat requiremets and temp leaving cooling pads (W)
    CoolE = VentRateReq/ventEff * 3600/1000; %energy requirements of ventilation system (kW)   
end

%function for finding electricity requirements of cooling system
%inputs: supplemental energy required, relative humidity, Temp outside,
%effficiency of evaporative cooling pad, heat capacity of air, density of
%air, minimum ventilation rate, set temperature, efficiency of ventilation
%system
%output: electricity requirements in kW*hr
function evapCoolE = evapCool(qSup, RH, TdryIn, evapPadEff, CpAir, rhoAir, ventRate, Tset, ventEff)
    TwetIn = TdryIn*atan(0.151977*((RH+8.313659)^0.5)) + atan(TdryIn+RH) - atan(RH - 1.676331) + 0.00391838*(RH^(3/2))*atan(0.023101*RH) - 4.686035;
    %^Calculation of wet bulb temperature
    TdryOut = TdryIn - evapPadEff*(TdryIn - TwetIn); %calculation of temp of air leaving evap cooling pads and going into the house
    qSupNoVent = qSup + ventE(CpAir, rhoAir, ventRate, Tset, TdryIn); %removal of ventilation from supplemental heating required because it is about to be recalculated (W)
    VentRateReq = qSupNoVent/(CpAir*rhoAir*(Tset - TdryOut)); %Ventilation rate required given supplemental heat requiremets and temp leaving cooling pads (W)
    evapCoolE = VentRateReq/ventEff * 3600/1000; %energy requirements of ventilation system (kW)   
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
function ventEnergy = ventE(CpAir, rhoAir, ventRate, Tset, Tout)
    ventEnergy = CpAir*rhoAir*ventRate*(Tset - Tout); % energy loss to ventilation
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
    for i = 1:1
        jan = sum(enrgCost(i:end,1:(24*31)), 'all');
        feb = sum(enrgCost(i:end,(24*31+1):(24*59)), 'all');
        mar = sum(enrgCost(i:end,(24*59+1):(24*90)), 'all');
        apr = sum(enrgCost(i:end,(24*90+1):(24*120)), 'all');
        may = sum(enrgCost(i:end,(24*120+1):(24*151)), 'all');
        jun = sum(enrgCost(i:end,(24*151+1):(24*181)), 'all');
        jul = sum(enrgCost(i:end,(24*181+1):(24*212)), 'all');
        aug = sum(enrgCost(i:end,(24*212+1):(24*243)), 'all');
        sep = sum(enrgCost(i:end,(24*243+1):(24*273)), 'all');
        oct = sum(enrgCost(i:end,(24*273+1):(24*304)), 'all');
        nov = sum(enrgCost(i:end,(24*304+1):(24*334)), 'all');
        dec = sum(enrgCost(i:end,(24*334+1):(24*365)), 'all');
%         fprintf('The average energy usage in January was %d kWh.\n', jan, 'all')/4)
%         fprintf('The average energy usage in February was %d kWh.\n', feb, 'all')/4)
%         fprintf('The average energy usage in March was %d kWh.\n', mar, 'all')/4)
%         fprintf('The average energy usage in April was %d kWh.\n', apr, 'all')/4)
%         fprintf('The average energy usage in May was %d kWh.\n', may, 'all')/4)
%         fprintf('The average energy usage in June was %d kWh.\n', jun, 'all')/4)
%         fprintf('The average energy usage in July was %d kWh.\n', jul, 'all')/4)
%         fprintf('The average energy usage in August was %d kWh.\n', aug, 'all')/4)
%         fprintf('The average energy usage in September was %d kWh.\n', sep, 'all')/4)
%         fprintf('The average energy usage in October was %d kWh.\n', oct, 'all')/4)
%         fprintf('The average energy usage in November was %d kWh.\n', nov, 'all')/4)
%         fprintf('The average energy usage in December was %d kWh.\n', dec, 'all')/4)
        figure(i+8)
        bar([jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec])
        title('Monthly Energy Usage for 2019')
        xlabel('Month')
        ylabel('Energy Usage (kWh)')
        set(gca,'xticklabel',{'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'})
    end
end
    
function yearlyUse(enrgCost)
    enrgCost = enrgCost / 1000; % convert to kWh
    perYearUse = sum(enrgCost, 2);
    fprintf('The energy usage in 2016 was %d kWh.\n', perYearUse(1));
    fprintf('The energy usage in 2017 was %d kWh.\n', perYearUse(2));
    fprintf('The energy usage in 2018 was %d kWh.\n', perYearUse(3));
    fprintf('The energy usage in 2019 was %d kWh.\n', perYearUse(4));
end
