%Get parameter info and vary output

function [FlywheeldiaOVary,torqueVary,powerVary] = getParamVary(Pmin,volumeC,volumeE,volumeR,Tc,theta2,length)
newTe = 1200;
[P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,newTe,theta2);
[Fp]  = getFp(P);
[torqueVary]  = getTorque(Fp,length,theta2 );
volumeT = volumeC+volumeE+volumeR;
plot(volumeT,P)
end