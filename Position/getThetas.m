function [ theta_0, theta_f, w_0]  = getThetas(T, T_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getThetas
%
%  PURPOSE
%   find the intersection of the engine torque and the average torque to
%   get theta_0, theta_f, and w_0
%  INPUTS
%   T: engine torque as a function of crank angle (N*m)
%   T_avg: average torque for one cycle (N*m)
%  OUTPUTS
%   theta_0: crank angle where energy is being added to the flywheel (deg)
%   theta_f: crank angle where energy is removed from the flywheel (deg)
%	w_0: angular velocity where energy is being added to the flywheel (rad)
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
theta_0 = fzero(fun,[0.1 180]);
theta_f = fzero(fun,[180 360]);


end

function [ T_diff ]  = torqueDiff(theta2, T, T_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: torqueDiff
%
%  PURPOSE
%   setup a function of the difference between the engine torque and
%   average torque
%  INPUTS
%   theta2: crank angle of the engine (deg)
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
%   index: location in the theta array corresponding to the input angle
%  FUNCTIONS CALLED
%   round: rounds to the nearest integer
%  START OF EXECUTABLE CODE
%
% find the index in the torque array corresponding to the input theta value
index = round(theta2/(360/3600));

% use the index above to calculate the difference in torques
T_diff = T(index)-T_avg;
end