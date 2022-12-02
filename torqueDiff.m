function [ T_diff ]  = torqueDiff(theta2, T, T_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: torqueDiff
%
%  PURPOSE
%   setup a function of the difference between the engine torque and
%   average torque
%  INPUTS
%   theta2: crank angle of the engine (rad)
%   T: engine torque as a function of theta (N*m)
%   T_avg: average torque over one cycle (N*m)
%  OUTPUT
%   T_diff: difference between engine torque and average torque to be set
%   to zero to find theta_0 and theta_f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/2/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%
%  FUNCTIONS CALLED
%
%  START OF EXECUTABLE CODE
%
% find the index in the torque array corresponding to the input theta value
index = round(theta2/(2*pi/3600));

% use the index above to calculate the difference in torques
T_diff = T(index)-T_avg;
end