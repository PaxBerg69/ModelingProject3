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
%	diffEQ: function handle defining the differential equation to be
%	used in ode45
%	ode45: differential equation solver
%	min: minimum value of an array
%	max: maximum value of an array
%
%  START OF EXECUTABLE CODE
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

% calculate the coefficient of fluctuation for our engine
COF_act = (w2_Max-w2_Min)/w_avg;
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