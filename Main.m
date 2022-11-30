%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: Project 3: Main Function
%  PURPOSE 
% 
%  INPUT
%  NA
%
%  OUTPUT
%  Print Output values
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
%Cylinder Values
cyl.stroke = 0.07; %Stroke [mm]
cyl.bore = 1;
cyl.CR = 1.58;

%Crank Values
crank.angle = 90;
[volume1, volume2] = cylVolumeFun(cyl,crank);
