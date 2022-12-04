function[difference] = flySize(I,Density,Width,ri)
% from the moment of inertia equation of a hollow cylinder
difference = 2*I/(Density*Width*pi) - 0.00002401 - 0.001372*ri - 0.0588*ri^2 - 0.28*ri^3;
end   % ends difference function