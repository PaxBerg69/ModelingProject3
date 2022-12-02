function [ theta_0, theta_f ]  = getThetas(T, T_avg, theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getThetas
%
%  PURPOSE
%   find the intersection of the engine torque and the average torque to
%   get theta_0 and theta_f
%  INPUTS
%   T: torque as a function of crank angle (N*m)
%   T_avg: average torque for one cycle (N*m)
%   theta: array of crank angles
%  OUTPUTS
%   theta_0: crank angle where energy is being added to the flywheel (rad)
%   theta_f: crank angle where energy is removed from the flywheel (rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/2/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   dif: function handle for the difference between engine torque and
%   average torque.
%  FUNCTIONS CALLED
%   
%  START OF EXECUTABLE CODE
% define the difference in engine torque and T_avg to use in fzero
% use a function handle for fzero

index_0 = fzero(dif,[0,theta(181)]);
index_f = fzero(dif,[theta(181),theta(361)]);

theta_0 = theta(index_0);
theta_f = theta(index_f);

end