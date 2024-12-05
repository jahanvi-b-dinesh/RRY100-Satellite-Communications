clc;

R_e = 6378.1;  %radius of earth
r = 42164.2; % semi-major axis
L = 57.7749; % latitude of ground station - Säve Airport (Säve Flygplats)
l = 303.7; % relative longitude of satellite 
h = R_e / r;

cos_phi = cosd(L) * cosd(l);
n = cos_phi * cos_phi;
disp(['cos_phi = ', num2str(cos_phi)]);

phi = acosd(cos_phi);
disp(['phi = ', num2str(phi)]);

%Distance of satellite from earth station
R = sqrt(R_e^2 + r^2 - (2*R_e*r*cos_phi)); 
disp(['R = ', num2str(R)]);

%Elevation angle
E = atand((cos_phi - h) / sqrt(1-n));
disp(['Elevation angle = ', num2str(E)]);

%Nadir Angle
Theta = asind(cosd(E) * h);
disp(['Nadir angle = ', num2str(Theta)]);

sum = phi + E + Theta;
disp(['sum = ',num2str(sum)]);