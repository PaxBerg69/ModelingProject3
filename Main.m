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

length.OaA = 0.0138;
length.OaC = 0.0138;
length.AB = 0.046;
length.CD = 0.0705;

theta3displacer = zeros(1,3600);
theta3power = zeros(1,3600);
theta2 = linspace(1,360,3600);
ydisplacer = zeros(1,3600);
ypower = zeros(1,3600);

% plot(theta2, theta3power, theta2, theta3displacer);
% legend('power piston', 'displacer piston');
% xlabel('theta2 (deg)');
% ylabel('theta3 (deg)');

%Crank Values
crank.angleP = 0 : 0.1 : 360;
crank.angleD = 90 : 0.1 : 450;
%Conrod Values
conRod.length = 0.055;

%% Function Calls
[theta3displacer, theta3power] = getTheta3(length,theta2);
[ydisplacer, ypower ]  = getYPosition( theta2, theta3displacer, theta3power, length);
[volumeE] = getVolumeE(cylD,ydisplacer);
% [volumeC] = getVolumeC()
plot(theta2,volumeE)
xlabel('Crank Angle [deg]')
ylabel('Volume [mm]')
hold
% plot(crank.angleD,volumeC)

