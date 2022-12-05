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
difference = 2*I/(Density*Width*pi) - 0.00002401 - 0.001372*ri - 0.0588*ri^2 - 0.28*ri^3;
end   % ends difference function