
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
