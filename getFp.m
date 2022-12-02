%DELETE COMMENTED PORTION BELOW IF PRESSURE IS A 3600 ITEM ARRAY AND AN FP
%ARRAY IS DEFINED IN MAIN


% pressure = linspace(0,200,3600);
% Fp = zeros(1,3600);
% 
% Fp = getFptest(pressure);

function [ Fp ]  = getFptest( pressure )
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
%
%  FUNCTIONS CALLED
%   None
%
%  START OF EXECUTABLE CODE
%

diam = 0.07;            %diameter of bore in meters
Fp = pressure * (3.14/4) * diam.^2;    %force on piston in N

end
