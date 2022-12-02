function [ V_c ]  = getVolumeC(disp, power, D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getVolumeC
%
%  PURPOSE
%   calculate the compression volume of the stirling engine
%  INPUTS
%   disp: structure containing information on the displacer piston
%   power: structure containing information on the power piston
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
V_c = 0.25*pi*(disp.y - power.y)*D^2;
end