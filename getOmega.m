function [w_2] = getOmega(Tavg,I_flywheel,T,Theta_0)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: Andrew Casar, Trey Weber
%  DATE:12/05/2022
%
%  DESCRIPTION OF LOCAL VARIABLES:
% 
%  FUNCTIONS CALLED: getdiffEQ 
%
%  START OF EXECUTABLE CODE
COF = 0.002;
w_avg = 2000;
fun = @(Theta_2,w_2) diffEQ(T,Tavg,I_flywheel,w_2,Theta_2);
w_0 = 2000*(2-COF)/2;
[Theta_2,w_2] = ode45(fun,[Theta_0 Theta_0+360],w_0);
w2_Min = min(w_2);       % in rad/s
w2_Max = max(w_2);       % in rad/s
if w2_Max <= COF*w_avg+w2_Min && w2_Min >= -COF*w_avg+w2_Max
    disp('The angular velocity of the crank is within the bounds')
else
    disp('The angular velocity of the crank is not within the bounds')
end  % ends if statement
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
%  DESCRIPTION OF LOCAL VARIABLES:
% 
%  FUNCTIONS CALLED: none
%
%  START OF EXECUTABLE CODE
% find the index in the torque array corresponding to the input theta value
index = round(Theta_2/(360/3600));

if index <= 3600
	dwdtheta = (T(index)-T_avg)/(I*w);   % the diffeq to find w
else
	dwdtheta = (T(index-3600)-T_avg)/(I*w);   % the diffeq to find w
end
end