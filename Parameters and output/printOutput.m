function [] = printOutput(torque,power,FlywheeldiaO,P,volumeT)
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

% Ideal
plot(volumeT,P);
xlabel('Volume [m]')
ylabel('Pressure [Pa]')

% % Actual
% plot(volumeT,Pactual);
% xlabel('Volume [m]')
% ylabel('Pressure [Pa]')

%Plot Torque, Power, Flywheel Diameter based on the varied parameter
% plot(torque,power,flywheeldiam,variedParam);
end