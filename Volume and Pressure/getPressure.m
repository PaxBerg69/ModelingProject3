function[P,Mtot] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getPressure
%
%  PURPOSE
%	calculate the pressure as a function of crank angle
%  INPUTS
%	Pmin: minimum pressure in the cylinder (at BDC) (Pa)
%	volumeC: compression volume as a function of crank angle (m^3)
%	volumeE: expansion volume as a function of crank angle (m^3)
%	volumeR: regenerator volume (m^3)
%	Tc: temperature in the compression region (K)
%	Te: temperature in the expansion region (K)
%	theta2: crank angle of the flywheel (deg)
%  OUTPUTS
%	P: pressure as a function of crank angle (Pa)
%	Mtot: total mass in the cylinder (kg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Andrew Casar
%  DATE: 12/7/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%	Tr: temperature in the regenerator (K)
%	R: ideal gas constant for air (Pa*m^3/kg*K)
%  FUNCTIONS CALLED
%	zeros: create an array of all zeros
%	length: return the length of an array
%  START OF EXECUTABLE CODE
% calculate the regenerator temperature
Tr =(Tc+Te)/2;
R = 287.039;   % Gas constant in PaM^3/KgK 

% Start at BDC, since we know pressure there and can calc Mtot
Mtot = (Pmin/R)*((volumeC(1)/Tc)+(volumeE(1)/Te)+(volumeR(1)/Tr));
% use the fact that Mtot doesn't change to find P at all crank angles
P = zeros(1,length(theta2));
for i = 1:length(theta2) % get Pressure for all values of Theta 2. 
    P(i) = (Mtot*R)/((volumeC(i)/Tc+volumeE(i)/Te+volumeR/Tr));
end   % Ends for loop
end   % Ends getPressure function 
