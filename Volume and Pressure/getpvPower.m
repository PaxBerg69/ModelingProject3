function [ pvPower, cycPower ]  = getpvPower( pressure, volumeT, Ptop, Pbottom, omega_avg )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: pvPower
%
%  PURPOSE: Extracts the power output from the P-V plot and compares it to
%  the idealized processes described by the stirling cycle
%
%  INPUT: pressure, volumeT, Ptop, Pbottom, omega_avg
%
%  OUTPUT: pvPower - power output from the p-v response in kW
%          cycPower - power output from the stirling cycle in kW
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/7/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%
%  FUNCTIONS CALLED
%
%  START OF EXECUTABLE CODE
%

% First determine the time per cycle of the engine
timePerCycle = (2*3.1415) / omega_avg;

% Work/time = power for the engine
workBounds = volumeT(1:3600);

workPV = trapz(workBounds, pressure);
pvPower = workPV / (1000 *timePerCycle);

% Work/time = power for the cycle - redefine bounds for the volume
newRange = linspace(min(volumeT), max(volumeT), 3600);

idealPV = trapz(newRange, Ptop) - trapz(newRange, Pbottom);
cycPower = idealPV / (1000 * timePerCycle);

end

