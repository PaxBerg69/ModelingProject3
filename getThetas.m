function [ theta_0, theta_f ]  = getThetas(T, T_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getThetas
%
%  PURPOSE
%   find the intersection of the engine torque and the average torque to
%   get theta_0 and theta_f
%  INPUTS
%   T: engine torque as a function of crank angle (N*m)
%   T_avg: average torque for one cycle (N*m)
%  OUTPUTS
%   theta_0: crank angle where energy is being added to the flywheel (rad)
%   theta_f: crank angle where energy is removed from the flywheel (rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/2/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   fun: function handle for the difference between engine torque and
%   average torque as a function of theta
%  FUNCTIONS CALLED
%   torqueDiff: function defining the difference in engine torque and
%   average torque
%   fzero: find the zero of a nonlinear function
%  START OF EXECUTABLE CODE
% define function handle to use in fzero
fun = @(theta2) torqueDiff(theta2,T,T_avg);

% use fzero to find the angle where the engine torque is equal to the
% average torque
theta_0 = fzero(fun,[0,pi]);
theta_f = fzero(fun,[pi,2*pi]);
end