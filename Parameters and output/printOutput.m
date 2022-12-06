function [] = printOutput(theta2,theta0,torque,power,FlywheeldiaO,P,volumeT,w_2, Pbot, Ptop, P1, P2, P3, P4)
% fprintf(torque);
% fprintf(power);
% fprintf(FlywheeldiaO);

%% Plotting
% plot(theta2, theta3power, theta2, theta3displacer);
% legend('power piston', 'displacer piston');
% xlabel('theta2 (deg)');
% ylabel('theta3 (deg)');

% plot(theta2,volumeE)
% xlabel('Crank Angle [deg]')
% ylabel('Volume [m]')
% hold
% plot(theta2,volumeC);
% legend('VolumeE', 'VolumeC');

minVol = min(volumeT);
maxVol = max(volumeT);
volRange = linspace(minVol, maxVol, 3600);

% Ideal
plot(volumeT,P);
hold on;
plot(volRange, Ptop, 'color','red')
plot(volRange, Pbot,'color','red');
xlabel('Volume [m]')
ylabel('Pressure [Pa]')
line([minVol, minVol],[P2,P1],'Color','R');
line([maxVol,maxVol],[P4,P3],'Color','R');
legend('Actual','Ideal');
hold off;

% % Actual
% plot(volumeT,Pactual);
% xlabel('Volume [m]')
% ylabel('Pressure [Pa]')

% plot w_2 vs theta2
figure
plot(theta2,w_2);
xlabel('Crank Angle (deg)');
ylabel('Angular Velocity of Flywheel (rad/s)');

%Plot Torque, Power, Flywheel Diameter based on the varied parameter
% plot(torque,power,flywheeldiam,variedParam);
end