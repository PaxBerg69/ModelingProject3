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

clc;
close all;

%% Initial Values
%Cylinder Values
cylP.stroke = 0.05; %Stroke [m]
cylP.bore = 0.07;
cylP.CR = 1.58;
cylD.stroke = 0.07;
cylD.bore = cylP.bore;
cylD.CR = cylP.CR;

%Link Values
length.OaA = 0.0138;
length.OaC = 0.0138;
length.AB = 0.046;
length.CD = 0.0705;

%Array Setup
theta3displacer = zeros(1,3600);
theta3power = zeros(1,3600);
theta2 = linspace(0.1,360,3600);
ydisplacer = zeros(1,3600);
ypower = zeros(1,3600);
Fp = zeros(1,3600);
torque = zeros(1,3600);
Pbot = zeros(1,3600);
Ptop = zeros(1,3600);

%Pressure values
Pmin = 500000; %[Pa] 
Tc = 300;
Te = 900;
volumeR = 0.00001; %[m]
theta0 = 0;
thetaF = 0;

% Flywheel values
Density = 1;
Width = 1;
ri = 1;
I = 1;
Cf = 0.002;
omega_avg = 2000*0.10472;

%Crank plotting Values
crank.angleP = 0 : 0.1 : 360;
crank.angleD = 90 : 0.1 : 450;

variedParam = 'input';

%% Function Calls
[theta3displacer, theta3power] = getTheta3(length,theta2);
[ydisplacer, ypower ]  = getYPosition( theta2, theta3displacer, theta3power, length);
[volumeE] = getVolumeE(cylD,ydisplacer);
[volumeC] = getVolumeC(ydisplacer, ypower, cylD);
volumeT = volumeE+volumeC;
[P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2);
[ P1, P2, P3, P4, Ptop, Pbot ]  = getIdeal( P, volumeC, volumeR, volumeE, Pbot, Ptop, Pmin );
Fp = getFp(P);
[torque] = getTorque(Fp,length,theta2);
Tavg = getTavg(theta2, torque);
[theta0, thetaF] = getThetas(torque, Tavg);
deltaKE = getDeltaKE(theta0, thetaF, Tavg, torque, theta2);
I = getI(deltaKE,Cf,omega_avg);
[flywheelDiaO] = getFlywheelsize(I);
w_2 = getOmega(Tavg,I,torque,theta0);
%[FlywheeldiaOVary,torqueVary,powerVary] = getParamVary(Pmin,volumeC,volumeE,volumeR,Tc,theta2,length);

printOutput(theta2,torque,power,flywheelDiaO,P,volumeT,w_2,Pbot,Ptop, P1, P2, P3, P4)

