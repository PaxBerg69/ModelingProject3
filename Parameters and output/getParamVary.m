%Get parameter info and vary output

function [FlywheeldiaO,torque,power] = getParamVary(Pmin,volumeT,Tc,theta2)
newTe = 400:1:1200;
[P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,newTe,theta2);
plot(volumeT,P)
end