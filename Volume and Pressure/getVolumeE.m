%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: Cylinder Volume Function
%  PURPOSE 
%  To calculate expansion volume based on cylinder geometry and crank
%  angle
%  INPUT
%  Cylinder geometry (CylD,ydisplacer), Crank angle
%
%  OUTPUT
%  Expansion volume as a functon of crank angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Paxton Berger
%  DATE: 11/30/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%  yEnd: Calculated distance from crank to cylinder head
%  FUNCTIONS CALLED
%  NA
%
%  START OF EXECUTABLE CODE
%
function [volumeE] = getVolumeE(cylD,ydisplacer)

yEnd = 0.1048;

volumeE = (yEnd-ydisplacer)*pi*0.25*cylD.bore^2;

end