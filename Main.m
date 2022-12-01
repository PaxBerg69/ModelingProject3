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
cylP.stroke = 0.05; %Stroke [m]
cylP.bore = 0.07;
cylP.CR = 1.58;
cylD.stroke = 0.07;
cylD.bore = cylP.bore;
cylD.CR = cylP.CR;

%Crank Values
crank.angleP = 0 : 0.1 : 360;
crank.angleD = 90 : 0.1 : 450;

%Conrod Values
conRod.length = 0.055;

[volumeP, volumeD] = CylinderVol(cylP,cylD,crank,conRod);

plot(crank.angleP,volumeP)
xlabel('Crank Angle [deg]')
ylabel('Volume [mm]')
hold
plot(crank.angleD,volumeD)

