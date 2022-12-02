function [ V_c ]  = getVolumeC(yDisp, yPower, D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getVolumeC
%
%  PURPOSE
%   calculate the compression volume of the stirling engine
%  INPUTS
%   yDisp: y position of the displacer piston
%   yPower: y position of the power piston (m)
%   D: cylinder diameter (m)
%  OUTPUT
%   V_c: compression volume of the stirling engine for all crank angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/2/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   none
%  FUNCTIONS CALLED
%   none
%  START OF EXECUTABLE CODE
%
V_c = 0.25*pi*(yDisp - yPower)*D^2;
end