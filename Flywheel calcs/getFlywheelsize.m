function[FlywheeldiaO] = getFlywheelsize(I)

density = 8000;  % in Kg/m^3.
width = 0.05;    % in meters

x0 = [0,10];  % want the positive value, 10m radius flywheel unlikely.
fun = @(ri) flySize(I,density,width,ri);
ri = fzero(fun,x0);  % in meters
ro=0.07+ri;    % in meters
mass = density*width*pi*(ro^2-ri^2); % mass of the flywheel
FlywheeldiaO = 2*ro; % Outer diameter of the flywheel
end    % Ends getFlywheelsize function
