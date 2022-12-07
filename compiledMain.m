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
Cf = 0.002;
omega_avg = 2000*0.10472;

%Crank plotting Values
crank.angleP = 0 : 0.1 : 360;
crank.angleD = 90 : 0.1 : 450;

%% Function Calls
[theta3displacer, theta3power] = getTheta3(length,theta2);
[ydisplacer, ypower ]  = getYPosition( theta2, theta3displacer, theta3power, length);
[volumeE] = getVolumeE(cylD,ydisplacer);
[volumeC] = getVolumeC(ydisplacer, ypower, cylD);
volumeT = volumeE+volumeC;
[P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2);
[P1,P2,P3,P4,Ptop,Pbot] = getIdeal(P,volumeC,volumeR,volumeE,Pbot,Ptop,Pmin);
Fp = getFp(P);
[torque] = getTorque(Fp,length,theta2, torque, theta3power);
Tavg = getTavg(theta2, torque);
[theta0, thetaF] = getThetas(torque, Tavg);
deltaKE = getDeltaKE(theta0, thetaF, Tavg, torque, theta2);
I = getI(deltaKE,Cf,omega_avg);
[flywheelDiaO] = getFlywheelsize(I);
[COF_act,w_2] = getOmega(Tavg,I,torque,theta2);
power = Tavg*omega_avg/1000; % in kW
printOutput(theta2,ydisplacer,ypower,torque,power,flywheelDiaO,P,volumeE,volumeC,volumeT,w_2,COF_act,Pbot,Ptop, P1, P2, P3, P4);

%Parameter Vary
figure;
hold on;
plotResolution = 25;
Thigh = linspace(400,2000,plotResolution);
powerV = zeros(1,plotResolution);
Dvary = zeros(1,plotResolution);
for j = 1:plotResolution
    Te = Thigh(j);
    [theta3displacer, theta3power] = getTheta3(length,theta2);
    [ydisplacer, ypower ]  = getYPosition( theta2, theta3displacer, theta3power, length);
    [volumeE] = getVolumeE(cylD,ydisplacer);
    [volumeC] = getVolumeC(ydisplacer, ypower, cylD);
    volumeT = volumeE+volumeC;
    [P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2);
    [P1,P2,P3,P4,Ptop,Pbot]  = getIdeal( P, volumeC, volumeR, volumeE, Pbot, Ptop, Pmin );
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

%% Function Definitions
function [theta3displacer, theta3power]  = getTheta3(length,theta2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getTheta3
%
%  PURPOSE: Return the angle of the connector link attached to both the
%  displacer the power piston as a function of theta2
%
%  INPUT: Relevant lengths, driver angle
%
%  OUTPUT: theta3displacer = angle displacer connector relative to GCS in degrees
%          theta3power = angle of power connector relative to GCS in
%          degrees
%          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/2/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%      vec - vector components (rotated or not denoted by star) of each link
%      theta3star - angle of link 3 (connector) relative to local
%                   coordinate frame
%
%  FUNCTIONS CALLED
%      None
%
%  START OF EXECUTABLE CODE
%
%% Establish original link length in vector form

thetaS = 90.0 * (6.28/360);

% For the displacer piston
for z = 1:3600

    vec.D = length.OaC * (cosd(theta2(z)) + i*sind(theta2(z)));   %define driving vector (link 2) in vector form

    vec.Dstar = vec.D * exp(-i*thetaS);
    vec.Dxstar = real(vec.Dstar);
    vec.Dystar = imag(vec.Dstar);
    vec.Cystar = -vec.Dystar;
    vec.Cxstar = (length.CD*length.CD - vec.Dystar*vec.Dystar)^0.5;   %define components of driver and connector vector through rotation of mechanism by thetaS

    theta3star = atand(vec.Cystar / vec.Cxstar);  
    theta3displacer(z) = theta3star + thetaS * (360/6.28);   %find angle of connector relative to GCS by rotating back by thetaS
end

%for the power piston
theta2disp = zeros(1,360);
for z = 1:3600
    theta2disp(z) = theta2(z) - 90;
    vec.D = length.OaA * (cosd(theta2disp(z)) + i*sind(theta2disp(z)));   %define driving vector (link 2) in vector form

    vec.Dstar = vec.D * exp(-i*thetaS);
    vec.Dxstar = real(vec.Dstar);
    vec.Dystar = imag(vec.Dstar);
    vec.Cystar = -vec.Dystar;
    vec.Cxstar = (length.AB*length.AB - vec.Dystar*vec.Dystar)^0.5;   %define components of driver and connector vector through rotation of mechanism by thetaS

    theta3star = atand(vec.Cystar / vec.Cxstar);  
    theta3power(z) = theta3star + thetaS * (360/6.28);   %find angle of connector relative to GCS by rotating back by thetaS
end


end
function [ydisplacer, ypower]  = getYPosition( theta2, theta3displacer, theta3power, length )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getYPosition
%
%  PURPOSE: Calculate the y position of the displacer and the power piston
%  relative to the ground pivot
%
%  INPUT: theta2 - angle of driver link in degrees
%         theta3displacer,power - angle of connector links in degrees
%         lengths of links in meters
%         
%
%  OUTPUT: ydisplacer - y position of displacer in meters
%          ypower - y position of power piston in meters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/2/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%
%  FUNCTIONS CALLED
%   None
%
%  START OF EXECUTABLE CODE
%  Determine y position of displacer and power piston:

for z = 1:3600 
    ydisplacer(z) = length.OaC * sind(theta2(z)) + length.CD * sind(theta3displacer(z));
    ypower(z) = length.OaA * sind(theta2(z)-90.0) + length.AB * sind(theta3power(z));

end

end
function [volumeE] = getVolumeE(cylD,ydisplacer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: Cylinder Volume Function
%  PURPOSE 
%  To calculate expansion volume based on cylinder geometry and crank
%  angle
%  INPUT
%  Cylinder geometry (CylD,ydisplacer), Crank angle
%
%  OUTPUT
%  Expansion volume as a functon of crank angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Paxton Berger
%  DATE: 11/30/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%  yEnd: Calculated distance from crank to cylinder head
%  FUNCTIONS CALLED
%  NA
%
%  START OF EXECUTABLE CODE
%

yEnd = 0.1048;

volumeE = (yEnd-ydisplacer)*pi*0.25*cylD.bore^2;

% cylC.Vd = cylC.stroke*pi*(cylC.bore^2)/4;
% cylC.Vc = (cylC.CR-1)/cylC.Vd;
% cylC.R = conRod.length/(cylC.stroke*0.5);
% 
% cylD.Vd = cylD.stroke*pi*(cylD.bore^2)/4;
% cylD.Vc = (cylD.CR-1)/cylD.Vd;
% cylD.R = conRod.length/(cylD.stroke*0.5);
% 
% volumeD = cylD.Vc*(1+(0.5)*(cylD.CR-1)*(cylD.R+1-cosd(crank.angleD)-((cylD.R^2)-sind(crank.angleD)).^0.5));
% volumeP = (cylD.stroke*pi*0.25*cylD.bore^2)-volumeD;
end
function [volumeC]  = getVolumeC1(ydisplacer, ypower, cylD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getVolumeC
%
%  PURPOSE
%   calculate the compression volume of the stirling engine
%  INPUTS
%   disp: structure containing information on the displacer piston
%   power: structure containing information on the power piston
%   D: cylinder diameter (m)
%  OUTPUT
%   V_c: compression volume of the stirling engine for all crank angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/2/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   none
%  FUNCTIONS CALLED
%   none
%  START OF EXECUTABLE CODE
%
volumeC = 0.25*pi*(ydisplacer - ypower)*cylD.bore^2;
end
function [P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2)
Tr =(Tc+Te)/2;
R = 287.039;   % Gas constant in PaM^3/KgK 
% Assume we start at BDC, since we know pressure there.
Mtot = (Pmin/R)*(volumeC(1)/Tc+volumeE(1)/Te+volumeR(1)/Tr);
P = zeros(1,length(theta2));
for i = 1:length(theta2) % get Pressure for all values of Theta 2. 
    P(i) = (Mtot*R)/((volumeC(i)/Tc+volumeE(i)/Te+volumeR/Tr));
end   % Ends for loop
end   % Ends getPressure function 
function [P1,P2,P3,P4,Ptop,Pbot] = getIdeal(volumeC,volumeR,volumeE,Pmin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getIdeal
%
%  PURPOSE: Calculate key points of the ideal plot as well as lines across
%  the top and bottom for high and low temperature
%
%  INPUT: pressure, total volume
%
%  OUTPUT: P1 - low volume high temperature pressure
%          P2 - low volume low temperature pressure
%          P3 - high volume high temperature pressure
%          P4 - high volume low temperature pressure  - all in Pa
%          Ptop, Pbottom - low and high temperature ideal pressures across
%                          each total volume value (in Pa)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/5/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%   R - gas constant for the ideal gas law
%   Thigh - high temperature of the cycle (Te) in K
%   Tlow - low temperature of the cycle (Tc) in K
%   Tr - temperature of the regenerator in K
%   Mtot - total mass in the system in kg
%   volumeT - total volume in the system in m^3
%   Vhigh, Vlow - max and min total volumes of the engine in m^3
%   volRange - array of volumes uniformly spaced between Vhigh and Vlow
%
%  FUNCTIONS CALLED
%
%  START OF EXECUTABLE CODE
%

%define total mass, temperatures, and R
R = 287.039;
Thigh = 900;
Tlow = 300;
Tr =(Tlow+Thigh)/2;
Mtot = (Pmin/R)*(volumeC(1)/Tlow+volumeE(1)/Thigh+volumeR(1)/Tr);

%calculate the total volume as well as the max and min
volumeT = volumeC + volumeE + volumeR;
Vhigh = max(volumeT);
Vlow = min(volumeT);

%find each corner of the ideal plot given combinations of high-low states
%of pressure and volume
P1 = Mtot * R * Thigh / Vlow;
P2 = Mtot * R * Tlow / Vlow;
P3 = Mtot * R * Thigh / Vhigh;
P4 = Mtot * R * Tlow / Vhigh;

%define a new volume range (from min to max) and fill idealized pressure
%arrays for each one - this will return the sloped portions of the
%idealized stirling cycle response at high (top) and low (bot) temperatures
VolRange = linspace(Vlow, Vhigh, 3600);
Ptop = Mtot * R * Thigh ./ VolRange;
Pbot = Mtot * R * Tlow ./ VolRange;

end
function [torque]  = getTorque( Fp, length, theta2, torque, theta3power )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getTorque
%
%  PURPOSE: Return the torque on the driver link necessary to 
%
%  INPUT - load force on power piston in N
%          length of crank arm of power piston in meters
%
%  OUTPUT - torque on crank arm as a function of driver angle (Nm)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/2/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%
%  FUNCTIONS CALLED
%   None
%
%  START OF EXECUTABLE CODE
%

length.AG3 = length.AB / 2.0;
length.BG3 = length.AB / 2.0;

F12x = zeros(1,3600);
F12y = zeros(1,3600);
F23x = zeros(1,3600);
F23y = zeros(1,3600);
F34x = zeros(1,3600);
F34y = zeros(1,3600);

for z = 1:3600
    A = [1 0 -1 0 0 0 0; 0 1 0 -1 0 0 0; 0 0 length.OaA*sind(theta2(z)) -length.OaA*cosd(theta2(z)) 0 0 1; 0 0 1 0 -1 0 0; 0 0 0 1 0 -1 0; 0 0 -length.AG3*sind(theta3power(z)-180) length.AG3*cosd(theta3power(z)-180) length.BG3*sind(theta3power(z)) -length.BG3*cosd(theta3power(z)) 0; 0 0 0 0 0 1 0];  %matrix with coefficients in front of each force/torque
    C = [0; 0; 0; 0; 0; 0; Fp(z)];   %matrix relating forces to weights and inertial forces/torques

    B = inv(A)*C;   %from linear algrebra - B matrix will fill with relevant force values

    torque(z) = -B(7);    %establish each force/torque component according to matrix setup (see pdf for list of the 7 equations used)
end
end
function [ T_avg ]  = getTavg(theta2,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getTavg
%
%  PURPOSE
%   compute the average torque as a function of crank angle (theta)
%  INPUTS
%   theta2: array of crank angles (rad)
%   T: torque as a function of crank angle (N*m)
%  OUTPUT
%   T_avg: average torque for one cycle (N*m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/2/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   none
%  FUNCTIONS CALLED
%   trapz: trapezoidal numerical integration
%  START OF EXECUTABLE CODE
% use the trapz function to calculate the average torque over the input
% theta array
T_avg = trapz(theta2,T)/(360.0);
end
function [ theta_0, theta_f, w_0]  = getThetas(T, T_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getThetas
%
%  PURPOSE
%   find the intersection of the engine torque and the average torque to
%   get theta_0, theta_f, and w_0
%  INPUTS
%   T: engine torque as a function of crank angle (N*m)
%   T_avg: average torque for one cycle (N*m)
%  OUTPUTS
%   theta_0: crank angle where energy is being added to the flywheel (deg)
%   theta_f: crank angle where energy is removed from the flywheel (deg)
%	w_0: angular velocity where energy is being added to the flywheel (rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/2/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   fun: function handle for the difference between engine torque and
%   average torque as a function of theta
%  FUNCTIONS CALLED
%   torqueDiff: function defining the difference in engine torque and
%   average torque
%   fzero: find the zero of a nonlinear function
%  START OF EXECUTABLE CODE
% define function handle to use in fzero
fun = @(theta2) torqueDiff(theta2,T,T_avg);

% use fzero to find the angle where the engine torque is equal to the
% average torque
theta_0 = fzero(fun,[0.1 180]);
theta_f = fzero(fun,[180 360]);


end
function [ T_diff ]  = torqueDiff(theta2, T, T_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: torqueDiff
%
%  PURPOSE
%   setup a function of the difference between the engine torque and
%   average torque
%  INPUTS
%   theta2: crank angle of the engine (deg)
%   T: engine torque as a function of theta (N*m)
%   T_avg: average torque over one cycle (N*m)
%  OUTPUT
%   T_diff: difference between engine torque and average torque to be set
%   to zero to find theta_0 and theta_f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/2/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   index: location in the theta array corresponding to the input angle
%  FUNCTIONS CALLED
%   round: rounds to the nearest integer
%  START OF EXECUTABLE CODE
%
% find the index in the torque array corresponding to the input theta value
index = round(theta2/(360/3600));

% use the index above to calculate the difference in torques
T_diff = T(index)-T_avg;
end
function [ deltaKE  ]  = getDeltaKE( theta0, thetaF, Tavg, torque, theta2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getDeltaKE
%
%  PURPOSE: Calculate the change in kinetic energy from the torque response
%  of the system
%
%  INPUT: torque, theta0, thetaF, theta2, Tavg
%
%  OUTPUT: deltaKE (J)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/5/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%
%  FUNCTIONS CALLED
%   trapz
%
%  START OF EXECUTABLE CODE
%

%convert to radians
theta2 = deg2rad(theta2);
theta0 = deg2rad(theta0);
thetaF = deg2rad(thetaF);

%set bounds of integration
theta0check = theta2 - theta0;
thetaFcheck = theta2 - thetaF;

[minimum,theta0index] = min(abs(theta0check));
[minimum,thetaFindex] = min(abs(thetaFcheck));

bounds = theta2(theta0index:thetaFindex);
deltaKE = trapz(bounds, torque(theta0index:thetaFindex)-Tavg);

end
function [ Iflywheel ]  = getI(KE, Cf, w_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getI
%
%  PURPOSE
%   calculate the moment of inertia of the flywheel of a stirling engine
%  INPUTS
%   KE: maximum change in kinetic energy of the flywheel (J)
%   Cf: coefficient of fluctuation of the flywheel angular velocity
%   w_avg: average angular velocity of the flywheel (rad/s)
%  OUTPUT
%   Iflywheel: mass moment of inertia of the flywheel (kg*m^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Trey Weber
%  DATE: 12/4/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%   none
%  FUNCTIONS CALLED
%   none
%  START OF EXECUTABLE CODE
% use the equation provided in lecture
Iflywheel = KE/(Cf*w_avg^2);
end
function[flywheelDiaO] = getFlywheelsize(I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getFlywheelsize
%
%  PURPOSE
%	calculate the flywheel dimensions for a given mass moment of inertia
%  INPUT
%	I: mass moment of inertia of the flywheel (kg*m^2)
%  OUTPUT
%	flywheelDiaO: outer diameter of the flywheel (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Andrew Casar
%  DATE: 12/5/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%	density: density of the flywheel (kg/m^3)
%	width: width of the flywheel (m)
%	x0: bounds of radius values to check in the fzero function call (m)
%	fun: function handle of the flywheel inertia equation
%	ri: internal radius of the flywheel (m)
%	ro: outer radius of the flywheel (m)
%  FUNCTIONS CALLED
%	Idif: function definition to be used in the fzero function call
%	fzero: calculate the zero of a nonlinear function
%  START OF EXECUTABLE CODE
%
density = 8000;  % in Kg/m^3.
width = 0.05;    % in meters

x0 = [0,1];  % want the positive value, 1m radius flywheel unlikely.
fun = @(ri) Idif(I,density,width,ri);
ri = fzero(fun,x0);  % in meters
ro=0.07+ri;    % in meters
flywheelDiaO = 2*ro; % Outer diameter of the flywheel
end    % Ends getFlywheelsize function
function[difference] = Idif(I,Density,Width,ri)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: Idif
%
%  PURPOSE
%	define the equation for the moment of inertia of a flywheel
%  INPUTS
%	I: mass moment of inertia of the flywheel (kg*m^2)
%	Density: density of the flywheel material (kg/m^3)
%	Width: width of the flywheel (m)
%	ri: internal radius of the flywheel (m)
%  OUTPUT
%	difference: placeholder value that will be minimized in the fzero
%	function call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Andrew Casar
%  DATE: 12/5/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%	none
%  FUNCTIONS CALLED
%	none
%  START OF EXECUTABLE CODE
% from the moment of inertia equation of a hollow cylinder
difference = (2*I/(Density*Width*pi)) - 0.00002401 - 0.001372*ri - 0.0588*ri^2 - 0.28*ri^3;
end   % ends difference function
function [COF_act,w_2] = getOmega(Tavg,I_flywheel,T,theta_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getCheck
%
%  PURPOSE: check that Coefficient of fluctuation is within 0.002
%
%  INPUTS:
%	Tavg: average torque for one cycle (N*m)
%	I_flywheel: mass moment of inertia of the flywheel (kg*m^2)
%	T: torque as a function of crank angle (N*m)
%	Theta_0: crank angle where energy is being added to the flywheel (deg)
%
%  OUTPUT:
%	w_2: angular velocity of the crank (rad/s)
%	COF_act: actual coefficient of fluctuation of our stirling engine
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Andrew Casar, Trey Weber
%  DATE:12/05/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%	COF: coefficient of fluctuation allowed for the flywheel
%	w_0: minimum angular velocity of the flywheel (rad/s)
%	theta_array: array from theta_0 to theta_0 on the next cycle used in
%				 ode45 (deg)
%	w2_min: minimum angular velocity of the flywheel (rad/s)
%	w2_max: maximum angular velocity of the flywheel (rad/s)
% 
%  FUNCTIONS CALLED 
%	getdiffEQ: function handle defining the differential equation to be
%	used in ode45
%	ode45: differential equation solver
%	min: minimum value of an array
%	max: maximum value of an array
%
%  START OF EXECUTABLE CODE
COF = 0.002;
w_avg = 2000*2*pi/60;

% convert to radians
theta_2 = deg2rad(theta_2);

% define the differential equation for ode45
fun = @(Theta_2,w_2) diffEQ(T,Tavg,I_flywheel,w_2,Theta_2);

% run an initial guess with initial condition at w_avg
[theta_2,w_2] = ode45(fun,theta_2,w_avg);
w_2 = w_2.'; % transpose to make it a row vector

% calculate the actual average of the initial guess
w_2avg = trapz(theta_2,w_2)/(2*pi);
offset = w_2avg-w_avg;

% make a new guess but offset the IC by the difference in averages
[theta_2,w_2] = ode45(fun,theta_2,w_avg-offset);
w_2 = w_2.'; % transpose to make it a row vector

w2_Min = min(w_2);       % in rad/s
w2_Max = max(w_2);       % in rad/s

COF_act = (w2_Max-w2_Min)/w_avg;
end  % ends the getCheck function
function [dwdtheta] = diffEQ(T,T_avg,I,w,Theta_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getdiffEQ
%
%  PURPOSE: create the equation for ode45 to use
%
%  INPUT: T,T_avg,I,w,Theta  
%
%  OUTPUT: dydt, the differential eqution to be used by ode45
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  AUTHOR: Andrew Casar, Trey Weber
%  DATE:12/05/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%	index: location in the torque array corresponding to the input theta
%  FUNCTIONS CALLED
%	round: round to the nearest integer
%  START OF EXECUTABLE CODE
% find the index in the torque array corresponding to the input theta value
index = round(Theta_2/(2*pi/3600));

if index <= 3600
	dwdtheta = (T(index)-T_avg)/(I*w);   % the diffeq to find w
else
	dwdtheta = (T(index-3600)-T_avg)/(I*w);   % the diffeq to find w
end
end
function [] = printOutput(theta2,ydisplacer,ypower,torque,power,FlywheeldiaO,P,volumeE,volumeC,volumeT,w_2,COF_act, Pbot, Ptop, P1, P2, P3, P4)
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

fprintf('Average Engine Power: %f (kW)\n',power);
fprintf('Flywheel Outer Diameter: %f (m)\n',FlywheeldiaO);
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
hold off

minVol = min(volumeT);
maxVol = max(volumeT);
volRange = linspace(minVol, maxVol, 3600);

% P-V diagram
figure;
plot(volumeT,P/1000);
hold on;
plot(volRange, Ptop/1000, 'color','red')
plot(volRange, Pbot/1000,'color','red');
xlabel('Volume [m^3]');
ylabel('Pressure [kPa]');
line([minVol, minVol],[P2,P1]/1000,'Color','R');
line([maxVol,maxVol],[P4,P3]/1000,'Color','R');
legend('Actual','Ideal');
title('P-V Diagram of the Engine vs. Ideal');
hold off;

% plot torque vs theta2
figure;
plot(theta2,torque);
xlabel('Crank Angle (deg)');
ylabel('Engine Torque (N*m)');
title('Engine Torque vs Crank Angle');

% plot w_2 vs theta2
figure;
hold on;
plot(theta2,w_2);
xlabel('Crank Angle (deg)');
ylabel('Angular Velocity of Flywheel (rad/s)');
title('Flywheel Angular Velocity vs Crank Angle');
ylim([208.5 210.5]);
xlim([0 360]);
end



