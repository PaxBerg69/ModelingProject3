
function [ P1, P2, P3, P4, Ptop, Pbot ]  = getIdeal( pressure, volumeC, volumeR, volumeE, Pbot, Ptop, Pmin )
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
%          P4 - high volume low temperature pressure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/5/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%
%  FUNCTIONS CALLED
%
%  START OF EXECUTABLE CODE
%

%define total mass and R
R = 287.039;
Thigh = 900;
Tlow = 300;
Tr =(Tlow+Thigh)/2;
Mtot = (Pmin/R)*(volumeC(1)/Tlow+volumeE(1)/Thigh+volumeR(1)/Tr);
volumeT = volumeC + volumeE + volumeR;
Vhigh = max(volumeT);
Vlow = min(volumeT);

P1 = Mtot * R * Thigh / Vlow;
P2 = Mtot * R * Tlow / Vlow;
P3 = Mtot * R * Thigh / Vhigh;
P4 = Mtot * R * Tlow / Vhigh;

VolRange = linspace(Vlow, Vhigh, 200);

Ptop = Mtot * R * Thigh ./ VolRange;
Pbot = Mtot * R * Tlow ./ VolRange;

end


