
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
%
%  FUNCTIONS CALLED
%   trapz
%
%  START OF EXECUTABLE CODE
%

%convert to radians
theta2 = theta2*(6.28/360);
theta0 = theta0*(6.28/360);
thetaF = thetaF*(6.28/360);

%set bounds of integration
theta0check = theta2 - theta0;
thetaFcheck = theta2 - thetaF;

[minimum,theta0index] = min(abs(theta0check));
[minimum,thetaFindex] = min(abs(thetaFcheck));

bounds = theta2(theta0index:thetaFindex);
deltaKE = trapz(bounds, torque(theta0index:thetaFindex)-Tavg);

end