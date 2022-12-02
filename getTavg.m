function [ T_avg ]  = getTavg(theta2,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getTavg
%
%  PURPOSE
%   compute the average torque as a function of crank angle (theta)
%  INPUTS
%   theta2: array of crank angles (rad)
%   T: torque as a function of crank angle (N*m)
%  OUTPUT
%   T_avg: average torque for one cycle (N*m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/2/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   none
%  FUNCTIONS CALLED
%   trapz: trapezoidal numerical integration
%  START OF EXECUTABLE CODE
% use the trapz function to calculate the average torque over the input
% theta array
T_avg = trapz(theta2,T)/(360.0);
end