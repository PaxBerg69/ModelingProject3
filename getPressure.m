function[P] = getPressure(Pmin,volumeC,volumeE,volumeR,Tc,Te,theta2)
Tr =(Tc+Te)/2;
R = 287.039;   % Gas constant in PaM^3/KgK 
% Assume we start at BDC, since we know pressure there.
Mtot = (Pmin/R)*(volumeC(1)/Tc+volumeE(1)/Te+volumeR(1)/Tr);
P = zeros(length(theta2));
for i = 1:length(theta2) % get Pressure for all values of Theta 2. 
    P(i) = (Mtot*R)/((volumeC(i)/Tc+volumeE(i)/Te+volumeR/Tr));
end   % Ends for loop
end   % Ends getPressure function 
