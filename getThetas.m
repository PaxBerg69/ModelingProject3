function [ theta_0, theta_f ]  = getThetas(theta2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getThetas
%
%  PURPOSE
%   find the intersection of the engine torque and the average torque to
%   get theta_0 and theta_f
%  INPUTS
%   theta2: array of crank angles (rad)
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
fun = @torqueDiff;

% use fzero to find the angle where the engine torque is equal to the
% average torque
theta_0 = fzero(fun,[0,theta2(1800)]);
theta_f = fzero(fun,[theta2(1800),theta2(3600)]);

end