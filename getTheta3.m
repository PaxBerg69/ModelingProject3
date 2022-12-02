%CAN DELETE COMMENTED PORTION BELOW IF VARIABLES ARE DEFINED IN MAIN


function [ theta3displacer, theta3power ]  = getTheta3(length,theta2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getTheta3
%
%  PURPOSE: Return the angle of the connector link attached to both the
%  displacer the power piston as a function of theta2
%
%  INPUT: Relevant lengths, driver angle
%
%  OUTPUT: theta3displacer = angle displacer connector relative to GCS in degrees
%          theta3power = angle of power connector relative to GCS in
%          degrees
%          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/2/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%      vec - vector components (rotated or not denoted by star) of each link
%      theta3star - angle of link 3 (connector) relative to local
%                   coordinate frame
%
%  FUNCTIONS CALLED
%      None
%
%  START OF EXECUTABLE CODE
%
%% Establish original link length in vector form

thetaS = 90.0 * (6.28/360);

% For the displacer piston
for z = 1:3600

    vec.D = length.OaC * (cosd(theta2(z)) + i*sind(theta2(z)));   %define driving vector (link 2) in vector form

    vec.Dstar = vec.D * exp(-i*thetaS);
    vec.Dxstar = real(vec.Dstar);
    vec.Dystar = imag(vec.Dstar);
    vec.Cystar = -vec.Dystar;
    vec.Cxstar = (length.CD*length.CD - vec.Dystar*vec.Dystar)^0.5;   %define components of driver and connector vector through rotation of mechanism by thetaS

    theta3star = atand(vec.Cystar / vec.Cxstar);  
    theta3displacer(z) = theta3star + thetaS * (360/6.28);   %find angle of connector relative to GCS by rotating back by thetaS
end

%for the power piston
theta2disp = zeros(1,360);
for z = 1:3600
    theta2disp(z) = theta2(z) + 90;
    vec.D = length.OaA * (cosd(theta2disp(z)) + i*sind(theta2disp(z)));   %define driving vector (link 2) in vector form

    vec.Dstar = vec.D * exp(-i*thetaS);
    vec.Dxstar = real(vec.Dstar);
    vec.Dystar = imag(vec.Dstar);
    vec.Cystar = -vec.Dystar;
    vec.Cxstar = (length.AB*length.AB - vec.Dystar*vec.Dystar)^0.5;   %define components of driver and connector vector through rotation of mechanism by thetaS

    theta3star = atand(vec.Cystar / vec.Cxstar);  
    theta3power(z) = theta3star + thetaS * (360/6.28);   %find angle of connector relative to GCS by rotating back by thetaS
end


end





