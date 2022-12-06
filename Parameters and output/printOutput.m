function [] = printOutput(theta2,ydisplacer,ypower,torque,power,FlywheeldiaO,P,volumeE,volumeC,volumeT,w_2, Pbot, Ptop, P1, P2, P3, P4)
%fprintf(torque);
%fprintf(power);
%fprintf(FlywheeldiaO);

%% Plotting
figure(1);
hold on
plot(theta2, ypower);
plot(theta2, ydisplacer);
legend('power piston', 'displacer piston');
xlabel('theta2 (deg)');
ylabel('theta3 (deg)');
xlim([0 360]);
title('Piston Positions vs Crank Angle');
hold off

figure(2);
hold on
plot(theta2,volumeE);
xlabel('Crank Angle [deg]')
ylabel('Volume [m^3]')
plot(theta2,volumeC);
plot(theta2,volumeT);
legend('VolumeE', 'VolumeC', 'Total Volume');
title('Region Volumes vs. Crank Angle');
xlim([0 360]);
hold off

minVol = min(volumeT);
maxVol = max(volumeT);
volRange = linspace(minVol, maxVol, 3600);

% Ideal
figure(3);
plot(volumeT,P);
hold on;
plot(volRange, Ptop, 'color','red')
plot(volRange, Pbot,'color','red');
xlabel('Volume [m^3]');
ylabel('Pressure [Pa]');
line([minVol, minVol],[P2,P1],'Color','R');
line([maxVol,maxVol],[P4,P3],'Color','R');
legend('Actual','Ideal');
title('P-V Diagram of the Engine vs. Ideal');
hold off;

% plot w_2 vs theta2
figure(4);
plot(theta2,w_2);
xlabel('Crank Angle (deg)');
ylabel('Angular Velocity of Flywheel (rad/s)');
title('Flywheel Angular Velocity vs Crank Angle');
xlim([0 360]);
end