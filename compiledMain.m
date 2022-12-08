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
%  Locals
%
%  OUTPUT
%  Average Engine Power (kW), Flywheel outer diameter (m), Engine
%  coefficient of fluctuation (-). Plot relevant parameters and trends
%  associated with increased hot temps.
%  AUTHOR: Paxton Berger, Andrew Casar, Carter Zehr, Trey Weber
%  DATE: 11/30/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%  Cylinder Strucutre (cylP,cylD): A set of values that describe the power
%  and displacer cylinder geometries
%  
%  Length Structure (OaA,OaC,AB,CD): A set of values that describe the
%  linkage geometries needed for thetas, y-postions, and volumes.
%
%  Pressure Values (Pmin,Tc,Te,volumeR): A set of values that describe the
%  internal ideal gas properties of the cylinders during motion
%
%  Angular Requirment values (Cf,omega_avg): A set of values that describe the
%  angular velocity parameters required for a properly functioning
%  flywheel.
%  
%  Crank Values: (crank.angleP, crank.angleD): Two arrays that represent
%  the possible theta values
%
%  FUNCTIONS CALLED
%  getTheta3, getYPosition, getVolumeE, getVolumeC, getPressure, getIdeal, getFp,
%  getTorque, getTavg, getThetas, torqueDiff, getDeltaKE, getI, Idif,
%  getOmega, diffEQ, getpvPower, printOutput, getParamVary.
%
%  START OF EXECUTABLE CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Angular Requirment values
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
% first try was barely within range, try a new value by increasing I
[COF_new,w_2new] = getOmega(Tavg,I+0.2,torque,theta2);
[flywheelDiaO_new] = getFlywheelsize(I+0.2);
power = Tavg*omega_avg/1000; % in kW
[pvPower, cycPower ] = getpvPower(P,volumeT,Ptop,Pbot,omega_avg );
printOutput(theta2,ydisplacer,ypower,torque,power,flywheelDiaO_new,P,Mtot,volumeE,volumeC,volumeT,w_2new,COF_new,Pbot,Ptop, P1, P2, P3, P4, pvPower, cycPower);
getParamVary(Pmin,volumeR,Tc,theta2,length,cylD,Cf,omega_avg,torque);

%% Function Definitions
function [theta3displacer, theta3power] = getTheta3(length,theta2)
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

thetaS = 90.0 * (6.28/360);   %angle of the S vector in degrees (pivot --> load)

% For the displacer piston - define driver (D) and connectors (C) in
% rotated domain
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

%for the power piston, do the same thing
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
function [ydisplacer, ypower] = getYPosition(theta2, theta3displacer, theta3power, length )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getYPosition
%
%  PURPOSE: Calculate the y position of the displacer and the power piston
%  relative to the ground pivot (meters)
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
%   None
%
%  FUNCTIONS CALLED
%   None
%
%  START OF EXECUTABLE CODE
%  Determine y position of displacer and power piston:

for z = 1:3600 
    ydisplacer(z) = length.OaC * sind(theta2(z)) + length.CD * sind(theta3displacer(z));
    ypower(z) = length.OaA * sind(theta2(z)-90.0) + length.AB * sind(theta3power(z));   %the y position can be determined by taking the y component of the driver and connector links and adding them together

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

end
function [volumeC] = getVolumeC(ydisplacer, ypower, cylD)
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
function [P,Mtot] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getPressure
%
%  PURPOSE
%	calculate the pressure as a function of crank angle
%  INPUTS
%	Pmin: minimum pressure in the cylinder (at BDC) (Pa)
%	volumeC: compression volume as a function of crank angle (m^3)
%	volumeE: expansion volume as a function of crank angle (m^3)
%	volumeR: regenerator volume (m^3)
%	Tc: temperature in the compression region (K)
%	Te: temperature in the expansion region (K)
%	theta2: crank angle of the flywheel (deg)
%  OUTPUTS
%	P: pressure as a function of crank angle (Pa)
%	Mtot: total mass in the cylinder (kg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Andrew Casar
%  DATE: 12/7/22
%
%  DESCRIPTION OF LOCAL VARIABLES
%	Tr: temperature in the regenerator (K)
%	R: ideal gas constant for air (Pa*m^3/kg*K)
%  FUNCTIONS CALLED
%	zeros: create an array of all zeros
%	length: return the length of an array
%  START OF EXECUTABLE CODE
% calculate the regenerator temperature
Tr =(Tc+Te)/2;
R = 287.039;   % Gas constant in PaM^3/KgK 

% Start at BDC, since we know pressure there and can calc Mtot
Mtot = (Pmin/R)*((volumeC(1)/Tc)+(volumeE(1)/Te)+(volumeR(1)/Tr));
% use the fact that Mtot doesn't change to find P at all crank angles
P = zeros(1,length(theta2));
for i = 1:length(theta2) % get Pressure for all values of Theta 2. 
    P(i) = (Mtot*R)/((volumeC(i)/Tc+volumeE(i)/Te+volumeR/Tr));
end   % Ends for loop
end   % Ends getPressure function 
function [P1,P2,P3,P4,Ptop,Pbot] = getIdeal(volumeC, volumeR, volumeE, Pmin )

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
function [Fp] = getFp(pressure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getFp
%
%  PURPOSE: Return force on the piston as a function of theta2 (N)
%
%  INPUT: pressure - pressure as a function of theta2 in Pa
%
%  OUTPUT: Fp - force on piston as a function of theta2 in N
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/2/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%     diam - diameter of bore in meters
%     Patm - atmospheric pressure in kPa
%
%  FUNCTIONS CALLED
%   None
%
%  START OF EXECUTABLE CODE
%

Patm = 101.3*1000;  %Patm in Pa
diam = 0.07;     %diameter of bore in meters
Fp = (pressure-Patm) * (pi/4) * diam.^2;    %force on piston in N is pressure times area (take difference between observed and actual

end
function [torque] = getTorque( Fp, length, theta2, torque, theta3power )
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
%   length.AG3, length.BG3 - distance from connector pins to center of mass
%   of the third link in meters
%   A, B, C - matrices for kinetostatic analysis
%
%  START OF EXECUTABLE CODE
%

length.AG3 = length.AB / 2.0;
length.BG3 = length.AB / 2.0;    %define distance from points A and B to the center of mass of the connector link (assumed to be the middle);

for z = 1:3600
    A = [1 0 -1 0 0 0 0; 0 1 0 -1 0 0 0; 0 0 length.OaA*sind(theta2(z)) -length.OaA*cosd(theta2(z)) 0 0 1; 0 0 1 0 -1 0 0; 0 0 0 1 0 -1 0; 0 0 -length.AG3*sind(theta3power(z)-180) length.AG3*cosd(theta3power(z)-180) length.BG3*sind(theta3power(z)) -length.BG3*cosd(theta3power(z)) 0; 0 0 0 0 0 1 0];  %matrix with coefficients in front of each force/torque
    C = [0; 0; 0; 0; 0; 0; Fp(z)];   %matrix relating forces to weights and inertial forces/torques

    B = inv(A)*C;   %from linear algrebra - B matrix will fill with relevant force values

    torque(z) = -B(7);    %establish each force/torque component according to matrix setup (see pdf for list of the 7 equations used)
end
end
function [T_avg] = getTavg(theta2,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getTavg
%
%  PURPOSE
%   compute the average torque as a function of crank angle (theta)
%  INPUTS
%   theta2: array of crank angles (deg)
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
% convert theta2 to radians
theta2 = deg2rad(theta2);

% use the trapz function to calculate the average torque over the input
% theta array
T_avg = trapz(theta2,T)/(2*pi);
end
function [theta_0, theta_f] = getThetas(T, T_avg)
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
function [T_diff] = torqueDiff(theta2, T, T_avg)
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
function [deltaKE] = getDeltaKE( theta0, thetaF, Tavg, torque, theta2)
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
%   theta2, theta0, thetaF - redefinitions of angles to radians
%   theta0check, thetaFcheck - array of thetas (theta2 - thetaX)
%   theta0index, thetaFindex - index of theta0 and thetaF in the theta2
%                              array
%   bounds - section of theta2 from theta0 to thetaF
%
%  FUNCTIONS CALLED
%   trapz
%
%  START OF EXECUTABLE CODE
%

%convert to radians
theta2 = deg2rad(theta2);
theta0 = deg2rad(theta0);
thetaF = deg2rad(thetaF);   %convert relevant angles to radians

%set bounds of integration
theta0check = theta2 - theta0;
thetaFcheck = theta2 - thetaF;

[minimum,theta0index] = min(abs(theta0check));
[minimum,thetaFindex] = min(abs(thetaFcheck));   %determine indices of theta0 and thetaF within the theta2 array

bounds = theta2(theta0index:thetaFindex);    %define subsection of theta2 corresponding to theta0-thetaF domain
deltaKE = trapz(bounds, torque(theta0index:thetaFindex)-Tavg);    %use trapz to determine area under the curve (kinetic energy absorbed by the flywheel)

end
function [Iflywheel] = getI(KE, Cf, w_avg)
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
difference = (2*I/(Density*Width*pi)) - 0.00002401 - 0.001372*ri - 0.0294*ri^2 - 0.28*ri^3;
end   % ends difference function
function [COF_act,w_2i] = getOmega(Tavg,I_flywheel,T,theta_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getCheck
%
%  PURPOSE: check that Coefficient of fluctuation is within 0.002
%
%  INPUTS:
%	Tavg: average torque for one cycle (N*m)
%	I_flywheel: mass moment of inertia of the flywheel (kg*m^2)
%	T: torque as a function of crank angle (N*m)
%	Theta_2: crank angle of the flywheel (deg)
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
%	w_avg: average angular velocity of the flywheel (rad/s)
%	w_2i: initial guess of w_2 using w_avg as the initial condition for
%	ode45 (rad/s)
%	w2i_avg: average value of the initial guess (rad/s)
%	w2_min: minimum angular velocity of the flywheel (rad/s)
%	w2_max: maximum angular velocity of the flywheel (rad/s)
%	w2_avg: average value of the angular velocity (rad/s)
%	offset: difference in the actual average to the specified average
%	angular velocity (rad/s)
%  FUNCTIONS CALLED 
%	diffEQ: function handle defining the differential equation to be
%	used in ode45
%	ode45: differential equation solver
%	min: minimum value of an array
%	max: maximum value of an array
%	trapz: trapezoidal numerical integration
%  START OF EXECUTABLE CODE
w_avg = 2000*2*pi/60;

% convert to radians
theta_2 = deg2rad(theta_2);

% define the differential equation for ode45
fun = @(Theta_2,w_2) diffEQ(T,Tavg,I_flywheel,w_2,Theta_2);

% run an initial guess with initial condition at w_avg
[theta_2,w_2i] = ode45(fun,theta_2,w_avg);
w_2i = w_2i.'; % transpose to make it a row vector

% calculate the actual average of the initial guess
w2i_avg = trapz(theta_2,w_2i)/(2*pi);
offset = w2i_avg-w_avg;

% make a new guess but offset the IC by the difference in averages
[theta_2,w_2] = ode45(fun,theta_2,w_avg-offset);
w_2 = w_2.'; % transpose to make it a row vector

% calculate the new average
w2_avg = trapz(theta_2,w_2)/(2*pi);

w2_Min = min(w_2);       % in rad/s
w2_Max = max(w_2);       % in rad/s

% calculate the coefficient of fluctuation for our engine
COF_act = (w2_Max-w2_Min)/w2_avg;
end  % ends the getCheck function
function [dwdtheta] = diffEQ(T,T_avg,I,w,Theta_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: diffEQ
%
%  PURPOSE: create the sum of torques equation for ode45 to use
%
%  INPUT: 
%	T: torque of the engine (N*m)
%	T_avg: average torque of the engine (N*m)
%	I: mass moment of inertia of the flywheel (kg*m^2)
%	w: rotational velocity of the flywheel (rad/s)
%	Theta_2: crank angle of the engine (rad)  
%
%  OUTPUT: 
%	dwdtheta: derivative of omega with respect to theta_2
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

% make sure that the index doesn't exceed the bounds
if index <= 3600
	dwdtheta = (T(index)-T_avg)/(I*w);   % the diffeq to find w
else
	dwdtheta = (T(index-3600)-T_avg)/(I*w);   % the diffeq to find w
end
end
function [ pvPower, cycPower ] = getpvPower( pressure, volumeT, Ptop, Pbottom, omega_avg )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: pvPower
%
%  PURPOSE: Extracts the power output from the P-V plot and compares it to
%  the idealized processes described by the stirling cycle
%
%  INPUT: pressure, volumeT, Ptop, Pbottom, omega_avg
%
%  OUTPUT: pvPower - power output from the p-v response in kW
%          cycPower - power output from the stirling cycle in kW
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Carter Zehr
%  DATE: 12/7/2022
%
%  DESCRIPTION OF LOCAL VARIABLES
%
%  FUNCTIONS CALLED
%
%  START OF EXECUTABLE CODE
%

% First determine the time per cycle of the engine
timePerCycle = (2*3.1415) / omega_avg;

% Work/time = power for the engine
workBounds = volumeT(1:3600);

workPV = trapz(workBounds, pressure);
pvPower = workPV / (1000 *timePerCycle);

% Work/time = power for the cycle - redefine bounds for the volume
newRange = linspace(min(volumeT), max(volumeT), 3600);

idealPV = trapz(newRange, Ptop) - trapz(newRange, Pbottom);
cycPower = idealPV / (1000 * timePerCycle);

end
function printOutput(theta2,ydisplacer,ypower,torque,power,FlywheeldiaO,P,Mtot,volumeE,volumeC,volumeT,w_2,COF_act, Pbot, Ptop, P1, P2, P3, P4, pvPower, cycPower)
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

fprintf('Flywheel Outer Diameter: %f m\n',FlywheeldiaO);
fprintf('Engine Coefficient of Fluctuation: %f\n',COF_act);
%% Plotting

% Piston Positions vs Crank Angle
figure;
hold on
plot(theta2, ypower);
plot(theta2, ydisplacer);
legend('power piston', 'displacer piston');
xlabel('Crank Angle (deg)');
ylabel('y Position (m)');
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

fprintf('\nPower from the average torque and angular speed: %f kW ',power);
fprintf('\nPower from the P-V plot: %f kW ',pvPower);
fprintf('\nIdealized power from the stirling cycle: %f kW ',cycPower);
fprintf('\nPV power as a percent of idealized power: %f percent\n',100*pvPower/cycPower);
end
function getParamVary(Pmin,volumeR,Tc,theta2,length,cylD,Cf,omega_avg,torque)
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
