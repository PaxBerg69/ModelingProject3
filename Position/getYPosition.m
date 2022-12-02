function [ ydisplacer, ypower ]  = getYPosition( theta2, theta3displacer, theta3power, length )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getYPosition
%
%  PURPOSE: Calculate the y position of the displacer and the power piston
%  relative to the ground pivot
%
%  INPUT: theta2 - angle of driver link in degrees
%         theta3displacer,power - angle of connector links in degrees
%         lengths of links in meters
%         
%
%  OUTPUT: ydisplacer - y position of displacer in meters
%          ypower - y position of power piston in meters
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
%  Determine y position of displacer and power piston:

for z = 1:3600 
    ydisplacer(z) = length.OaC * sind(theta2(z)) + length.CD * sind(theta3displacer(z));
    ypower(z) = length.OaA * sind(theta2(z)-90.0) + length.AB * sind(theta3power(z));

end

end










