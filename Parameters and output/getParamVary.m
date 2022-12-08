%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getParamVary
%  PURPOSE
%  To change a critical parameter of the engine and determine the effects
%  on the torque, power, and flywheel geometry.
% 
%  INPUT
%  All info similar to main output function, as well varied parameter
%  selection.
%
%  OUTPUT
%  Average Engine Power (kW), Flywheel outer diameter (m) plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Paxton Berger
%  DATE: 12/6/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%  plotResolution, a variable determining the number of plotted points to
%  show trend of varied parameter agaisnt output values.
%
%  FUNCTIONS CALLED
%  getTheta3,getYPosition,gotVolumeE,getVolumeC,getPressure,getIdeal,
%  getFp,getTorque,getTavg,getThetas,getDeltaKE,getI,getFlywheelsize.getOmega.
%
%  START OF EXECUTABLE CODE

function getParamVary(Pmin,volumeR,Tc,theta2,length,cylD,Cf,omega_avg)
%% Parameter Vary
figure;
hold on;
plotResolution = 25; %Number of points on plots
Thigh = linspace(400,2000,plotResolution); %Setup varied parameter array space.
powerV = zeros(1,plotResolution); %Initialize output arrays
Dvary = zeros(1,plotResolution); %Initialize output arrays
for j = 1:plotResolution %For loop to iterate varied parameter and store in array
    Te = Thigh(j);
    [theta3displacer, theta3power] = getTheta3(length,theta2);
    [ydisplacer, ypower ]  = getYPosition( theta2, theta3displacer, theta3power, length);
    [volumeE] = getVolumeE(cylD,ydisplacer);
    [volumeC] = getVolumeC(ydisplacer, ypower, cylD);
    [P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2);
    Fp = getFp(P);
    [torque] = getTorque(Fp,length,theta2, torque, theta3power);
    Tavg = getTavg(theta2, torque);
    [theta0, thetaF] = getThetas(torque, Tavg);
    deltaKE = getDeltaKE(theta0, thetaF, Tavg, torque, theta2);
    I = getI(deltaKE,Cf,omega_avg);
    Dvary(j) = 1000*getFlywheelsize(I);
    powerV(j) = Tavg*omega_avg/1000;
end
%Plot
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
end