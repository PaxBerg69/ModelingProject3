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
function [volumeP, volumeD] = getVolumeE(cylC,cylD,crank,conRod)
cylC.Vd = cylC.stroke*pi*(cylC.bore^2)/4;
cylC.Vc = (cylC.CR-1)/cylC.Vd;
cylC.R = conRod.length/(cylC.stroke*0.5);

cylD.Vd = cylD.stroke*pi*(cylD.bore^2)/4;
cylD.Vc = (cylD.CR-1)/cylD.Vd;
cylD.R = conRod.length/(cylD.stroke*0.5);

volumeD = cylD.Vc*(1+(0.5)*(cylD.CR-1)*(cylD.R+1-cosd(crank.angleD)-((cylD.R^2)-sind(crank.angleD)).^0.5));
volumeP = (cylD.stroke*pi*0.25*cylD.bore^2)-volumeD;




end