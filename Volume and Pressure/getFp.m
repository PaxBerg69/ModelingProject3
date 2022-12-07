function [ Fp ]  = getFp( pressure )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getFp
%
%  PURPOSE: Return force on the piston as a function of theta2 (N)
%
%  INPUT: pressure - pressure as a function of theta2 in Pa
%
%  OUTPUT: Fp - force on piston as a function of theta2 in N
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/2/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%     diam - diameter of bore in meters
%     Patm - atmospheric pressure in kPa
%
%  FUNCTIONS CALLED
%   None
%
%  START OF EXECUTABLE CODE
%

Patm = 101.3*1000;  %Patm in Pa
diam = 0.07;     %diameter of bore in meters
Fp = (pressure-Patm) * (pi/4) * diam.^2;    %force on piston in N is pressure times area (take difference between observed and actual

end
