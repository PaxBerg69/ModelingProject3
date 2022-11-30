%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: Cylinder Volume Function
%  PURPOSE 
%  To calculate both cylinder volumes based on cylinder geometry and crank
%  angle
%  INPUT
%  Crank angle
%
%  OUTPUT
%  Volume as a functon of crank angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Paxton Berger
%  DATE: 11/30/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%  
%  FUNCTIONS CALLED
%  NA
%
%  START OF EXECUTABLE CODE
%
function [volume1, volume2] = CylinderVol(cyl,crank)
Vd = cyl.stroke*pi*(cyl.bore^2)/4;
Vc = (cyl.CR-1)/Vd

end