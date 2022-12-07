function [ deltaKE  ]  = getDeltaKE( theta0, thetaF, Tavg, torque, theta2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getDeltaKE
%
%  PURPOSE: Calculate the change in kinetic energy from the torque response
%  of the system
%
%  INPUT: torque, theta0, thetaF, theta2, Tavg
%
%  OUTPUT: deltaKE (J)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/5/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%   theta2, theta0, thetaF - redefinitions of angles to radians
%   theta0check, thetaFcheck - array of thetas (theta2 - thetaX)
%   theta0index, thetaFindex - index of theta0 and thetaF in the theta2
%                              array
%   bounds - section of theta2 from theta0 to thetaF
%
%  FUNCTIONS CALLED
%   trapz
%
%  START OF EXECUTABLE CODE
%

%convert to radians
theta2 = deg2rad(theta2);
theta0 = deg2rad(theta0);
thetaF = deg2rad(thetaF);   %convert relevant angles to radians

%set bounds of integration
theta0check = theta2 - theta0;
thetaFcheck = theta2 - thetaF;

[minimum,theta0index] = min(abs(theta0check));
[minimum,thetaFindex] = min(abs(thetaFcheck));   %determine indices of theta0 and thetaF within the theta2 array

bounds = theta2(theta0index:thetaFindex);    %define subsection of theta2 corresponding to theta0-thetaF domain
deltaKE = trapz(bounds, torque(theta0index:thetaFindex)-Tavg);    %use trapz to determine area under the curve (kinetic energy absorbed by the flywheel)

end