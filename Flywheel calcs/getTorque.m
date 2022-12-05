% DELETE COMMENTED PORTION IF ALL THESE ARE DEFINED ELSEWHERE

% Fp = linspace(0,500,3600);
% torque = zeros(1,3600);
% theta2 = linspace(0.1,360,3600);
% length.OaA = 0.0138;
% torque = getTorquetest(Fp, length, theta2);
% 
% plot(theta2, torque);


function [torque,power]  = getTorque( Fp, length, theta2 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getTorque
%
%  PURPOSE: Return the torque on the driver link necessary to 
%
%  INPUT - load force on power piston in N
%          length of crank arm of power piston in meters
%
%  OUTPUT - torque on crank arm as a function of driver angle (Nm)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/2/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%
%  FUNCTIONS CALLED
%   None
%
%  START OF EXECUTABLE CODE
%

torque = -Fp .* cosd(theta2 - 90.0) .* length.OaA;
power = torque*2000;
end




