function [ Iflywheel ]  = getI(KE, Cf, w_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getI
%
%  PURPOSE
%   calculate the moment of inertia of the flywheel of a stirling engine
%  INPUTS
%   KE: maximum change in kinetic energy of the flywheel (J)
%   Cf: coefficient of fluctuation of the flywheel angular velocity
%   w_avg: average angular velocity of the flywheel (rad/s)
%  OUTPUT
%   Iflywheel: mass moment of inertia of the flywheel (kg*m^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/4/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   none
%  FUNCTIONS CALLED
%   none
%  START OF EXECUTABLE CODE
% use the equation provided in lecture
Iflywheel = KE/(Cf*w_avg^2);
end
