%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: Project 3: Main Function
%  PURPOSE
%  To model a Stirling Engine using the given design parameters through
%  volume, pressure, and rotational analysis, with a goal of deriving the
%  proper flywheel diameter required for operation, as well as the output
%  torque and power values, while maintaning an acceptable level of angular
%  velocity fluctuation.
% 
%  INPUT
%  NA
%
%  OUTPUT
%  Average Engine Power (kW), Flywheel outer diameter (m), Engine
%  coefficient of fluctuation (-). Plot relevant parameters and trends
%  associated with increased hot temps.
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
theta2 = linspace(0.1,360,3600);
torque = zeros(1,3600);

%Pressure values
Pmin = 500000; %[Pa] 
Tc = 300;
Te = 900;
volumeR = 0.00001; %[m]

% Flywheel values
Cf = 0.002;
omega_avg = 2000*2*pi/60;

%Crank plotting Values
crank.angleP = 0 : 0.1 : 360;
crank.angleD = 90 : 0.1 : 450;

%% Function Calls
[theta3displacer, theta3power] = getTheta3(length,theta2);
[ydisplacer, ypower ]  = getYPosition( theta2, theta3displacer, theta3power, length);
[volumeE] = getVolumeE(cylD,ydisplacer);
[volumeC] = getVolumeC(ydisplacer, ypower, cylD);
volumeT = volumeE+volumeC+volumeR;
[P,Mtot] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2);
[ P1, P2, P3, P4, Ptop, Pbot ]  = getIdeal(volumeC, volumeR, volumeE, Pmin );
Fp = getFp(P);
[torque] = getTorque(Fp,length,theta2, torque, theta3power);
Tavg = getTavg(theta2, torque);
[theta0, thetaF] = getThetas(torque, Tavg);
deltaKE = getDeltaKE(theta0, thetaF, Tavg, torque, theta2);
I = getI(deltaKE,Cf,omega_avg);
[flywheelDiaO] = getFlywheelsize(I);
[COF_act,w_2] = getOmega(Tavg,I,torque,theta2);
power = Tavg*omega_avg/1000; % in kW
[pvPower, cycPower ] = getpvPower( P, volumeT, Ptop, Pbot, omega_avg );
printOutput(theta2,ydisplacer,ypower,torque,power,flywheelDiaO,P,Mtot,volumeE,volumeC,volumeT,w_2,COF_act,Pbot,Ptop, P1, P2, P3, P4, pvPower, cycPower);

%% Parameter Vary
figure;
hold on;
plotResolution = 25;
Thigh = linspace(400,2000,plotResolution);
powerV = zeros(1,plotResolution);
Dvary = zeros(1,plotResolution);
for j = 1:plotResolution
<<<<<<< HEAD
    Te = Thigh(j);
    [theta3displacer, theta3power] = getTheta3(length,theta2);
    [ydisplacer, ypower ]  = getYPosition( theta2, theta3displacer, theta3power, length);
    [volumeE] = getVolumeE(cylD,ydisplacer);
    [volumeC] = getVolumeC(ydisplacer, ypower, cylD);
    volumeT = volumeE+volumeC;
    [P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2);
    [ P1, P2, P3, P4, Ptop, Pbot ]  = getIdeal( P, volumeC, volumeR, volumeE, Pbot, Ptop, Pmin );
    Fp = getFp(P);
    [torque] = getTorque(Fp,length,theta2, torque, theta3power);
    Tavg = getTavg(theta2, torque);
    [theta0, thetaF] = getThetas(torque, Tavg);
    deltaKE = getDeltaKE(theta0, thetaF, Tavg, torque, theta2);
    I = getI(deltaKE,Cf,omega_avg);
    Dvary(j) = 1000*getFlywheelsize(I);
    w_2 = getOmega(Tavg,I,torque,theta2);
    powerV(j) = Tavg*omega_avg/1000;
    [COF_act,w_2] = getOmega(Tavg,I,torque,theta2);
=======
Te = Thigh(j);
[theta3displacer, theta3power] = getTheta3(length,theta2);
[ydisplacer, ypower ]  = getYPosition( theta2, theta3displacer, theta3power, length);
[volumeE] = getVolumeE(cylD,ydisplacer);
[volumeC] = getVolumeC(ydisplacer, ypower, cylD);
volumeT = volumeE+volumeC;
[P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2);
[ P1, P2, P3, P4, Ptop, Pbot ]  = getIdeal( volumeC, volumeR, volumeE, Pmin );
Fp = getFp(P);
[torque] = getTorque(Fp,length,theta2, torque, theta3power);
Tavg = getTavg(theta2, torque);
[theta0, thetaF] = getThetas(torque, Tavg);
deltaKE = getDeltaKE(theta0, thetaF, Tavg, torque, theta2);
I = getI(deltaKE,Cf,omega_avg);
Dvary(j) = 1000*getFlywheelsize(I);
powerV(j) = Tavg*omega_avg/1000;
[COF_act,w_2] = getOmega(Tavg,I,torque,theta2);
>>>>>>> 696e184b0a932e4d82cef5a2adfd715b3f5cf95d
end
hold off

plot(Thigh, powerV);
title('Power Output as a Function of Te');
xlabel('Expansion Temperature (K)');
ylabel('Power Output (kW)');

figure;
plot(Thigh, Dvary);
title('Flywheel Diameter as a Function of Te');
xlabel('Expansion Temperature (K)');
ylabel('Flywheel Diameter (mm)');
