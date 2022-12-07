function [ P1, P2, P3, P4, Ptop, Pbot ]  = getIdeal( volumeC, volumeR, volumeE, Pmin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getIdeal
%
%  PURPOSE: Calculate key points of the ideal plot as well as lines across
%  the top and bottom for high and low temperature
%
%  INPUT: pressure, total volume
%
%  OUTPUT: P1 - low volume high temperature pressure
%          P2 - low volume low temperature pressure
%          P3 - high volume high temperature pressure
%          P4 - high volume low temperature pressure  - all in Pa
%          Ptop, Pbottom - low and high temperature ideal pressures across
%                          each total volume value (in Pa)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/5/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%   R - gas constant for the ideal gas law
%   Thigh - high temperature of the cycle (Te) in K
%   Tlow - low temperature of the cycle (Tc) in K
%   Tr - temperature of the regenerator in K
%   Mtot - total mass in the system in kg
%   volumeT - total volume in the system in m^3
%   Vhigh, Vlow - max and min total volumes of the engine in m^3
%   volRange - array of volumes uniformly spaced between Vhigh and Vlow
%
%  FUNCTIONS CALLED
%
%  START OF EXECUTABLE CODE
%

%define total mass, temperatures, and R
R = 287.039;
Thigh = 900;
Tlow = 300;
Tr =(Tlow+Thigh)/2;
Mtot = (Pmin/R)*(volumeC(1)/Tlow+volumeE(1)/Thigh+volumeR(1)/Tr);

%calculate the total volume as well as the max and min
volumeT = volumeC + volumeE + volumeR;
Vhigh = max(volumeT);
Vlow = min(volumeT);

%find each corner of the ideal plot given combinations of high-low states
%of pressure and volume
P1 = Mtot * R * Thigh / Vlow;
P2 = Mtot * R * Tlow / Vlow;
P3 = Mtot * R * Thigh / Vhigh;
P4 = Mtot * R * Tlow / Vhigh;

%define a new volume range (from min to max) and fill idealized pressure
%arrays for each one - this will return the sloped portions of the
%idealized stirling cycle response at high (top) and low (bot) temperatures
VolRange = linspace(Vlow, Vhigh, 3600);
Ptop = Mtot * R * Thigh ./ VolRange;
Pbot = Mtot * R * Tlow ./ VolRange;

end


