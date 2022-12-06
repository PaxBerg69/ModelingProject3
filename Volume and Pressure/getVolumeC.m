function [ volumeC ]  = getVolumeC(ydisplacer, ypower, cylD)
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
volumeC = 0.25*pi*(ydisplacer - ypower)*cylD.bore^2;
end
