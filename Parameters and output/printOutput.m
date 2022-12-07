%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: Print Output Function
%  PURPOSE 
%  Display are values and plots desired for analysis
%  INPUT
%  All relavent outputs (volume, torque, theta, power, location,
%  coefficient of fluctuation, varied parameters, etc.)
%
%  OUTPUT
%  No variable output,only print all values and plot relavent values/trends
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
function [] = printOutput(theta2,ydisplacer,ypower,torque,power,FlywheeldiaO,P,Mtot,volumeE,volumeC,volumeT,w_2,COF_act, Pbot, Ptop, P1, P2, P3, P4, pvPower, cycPower)
fprintf('Flywheel Outer Diameter: %f m\n',FlywheeldiaO);
fprintf('Engine Coefficient of Fluctuation: %f\n',COF_act);
%% Plotting

% Piston Positions vs Crank Angle
figure;
hold on
plot(theta2, ypower);
plot(theta2, ydisplacer);
legend('power piston', 'displacer piston');
xlabel('theta2 (deg)');
ylabel('theta3 (deg)');
xlim([0 360]);
title('Piston Positions vs Crank Angle');
hold off

% Volume vs Crank Angle
figure;
hold on
plot(theta2,volumeE);
xlabel('Crank Angle [deg]')
ylabel('Volume [m^3]')
plot(theta2,volumeC);
plot(theta2,volumeT);
legend('VolumeE', 'VolumeC', 'Total Volume');
title('Region Volumes vs. Crank Angle');
xlim([0 360]);
ylim([0 4e-04]);
hold off

minVol = min(volumeT);
maxVol = max(volumeT);
volRange = linspace(minVol, maxVol, 3600);

% P-V diagram (pressure converted to kPa, volume converted to specific)
figure;
plot(volumeT/Mtot,P/1000);
hold on;
plot(volRange/Mtot, Ptop/1000, 'color','red')
plot(volRange/Mtot, Pbot/1000,'color','red');
xlabel('Specific Volume [m^3/kg]');
ylabel('Pressure [kPa]');
line([minVol, minVol]/Mtot,[P2,P1]/1000,'Color','R');
line([maxVol,maxVol]/Mtot,[P4,P3]/1000,'Color','R');
legend('Actual','Ideal');
title('P-V Diagram of the Engine vs. Ideal');
xlim([0.14 0.28]);
ylim([200 1800]);
hold off;

% plot torque vs theta2
figure;
plot(theta2,torque);
xlabel('Crank Angle (deg)');
ylabel('Engine Torque (N*m)');
title('Engine Torque vs Crank Angle');
xlim([0 360]);

% plot w_2 vs theta2
figure;
hold on;
plot(theta2,w_2);
xlabel('Crank Angle (deg)');
ylabel('Angular Velocity of Flywheel (rad/s)');
title('Flywheel Angular Velocity vs Crank Angle');
yline(209.23);
yline(209.65);
ylim([208 211]);
xlim([0 360]);
legend('Flywheel Angular Velocity','Allowable Range');
hold off;

fprintf('\nPower from the average torque and angular speed - ');
fprintf('%f',power);
fprintf(' kW');
fprintf('\nPower from the P-V plot - ');
fprintf('%f',pvPower);
fprintf(' kW');
fprintf('\nIdealized power from the stirling cycle - ');
fprintf('%f',cycPower);
fprintf(' kW');
fprintf('\nPV power as a percent of idealized power - ');
fprintf('%f',100*pvPower/cycPower);
fprintf(' percent\n');
end