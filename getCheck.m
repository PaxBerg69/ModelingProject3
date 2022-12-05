function[] = getCheck(Tavg,I_flywheel,T,Theta_0,w_0,Theta_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: getCheck
%
%  PURPOSE: check that Coefficient of fluctuation is within 0.002
%
%  INPUT:Tavg,I_flywheel,T,Theta_0,w_0,Theta_2 
%
%  OUTPUT: 
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
fun = @(Theta_2,w_2) getdiffEQ(T,Tavg,I_flywheel,w_2,Theta_2);
[Theta_2,w_2] = ode45(fun,[Theta_0 Theta_0+360]);
w2_Min = min(w2);       % in rad/s
w2_Max = max(w2);       % in rad/s
w2_avg = (w2_Min+w2_Max)/2;  % in rad/s
if w2_Max <= COF*w2_avg+w2_Min && w2_Min >= -COF*w2_avg+w2_Max
    disp('The angular velocity of the crank is within the bounds')
else
    disp('The angular velocity of the crank is not within the bounds')
end  % ends if statement
end  % ends the getCheck function