clc;

RE = 6378; 
r = 42162.9;


l = 57.7749;
L1 = 5; %SAT1
L2 = 35; %B-SAT1
L3 = 57; %SAT2
L4 = 20; %B-SAT2


delta_theta_1 = L1 - L2;
delta_theta_2 = L2 - L3;
delta_theta_3 = L3 - L4;

% Calculate angular separation phi
cos_phi_1 = cosd(l) * cosd(L1);
phi_1 = acos(cos_phi_1);

cos_phi_2 = cosd(l) * cosd(L2);
phi_2 = acos(cos_phi_2); 

cos_phi_3 = cosd(l) * cosd(L3);
phi_3 = acos(cos_phi_3);

cos_phi_4 = cosd(l) * cosd(L4);
phi_4 = acosd(cos_phi_4); 

% Calculate elevation angle E
% numerator = cos_phi - RE / r;
% denominator = sqrt(1 - cos_phi^2);
% E = atan(numerator / denominator); % Elevation angle in radians


% E_deg = rad2deg(E);


%Distance from ground station to satellite

R1 = sqrt(RE^2 + r^2 - 2*RE*r*cosd(phi_1));
R2 = sqrt(RE^2 + r^2 - 2*RE*r*cosd(phi_2));
R3 = sqrt(RE^2 + r^2 - 2*RE*r*cosd(phi_3));
R4 = sqrt(RE^2 + r^2 - 2*RE*r*cosd(phi_4));
fprintf('Distance between ground station and Sat-1: %.2f kms\n', R1);
fprintf('Distance between ground station and B-Sat1: %.2f kms\n', R2);
fprintf('Distance between ground station and Sat-2: %.2f kms\n', R3);
fprintf('Distance between ground station and B-Sat2: %.2f kms\n', R4);


% Distance between satellites 
R12 = r * delta_theta_1;
R23 = r * delta_theta_2;
R34 = r * delta_theta_3;
fprintf('Distance between Sat-1 and B-Sat1: %.2f kms\n', R12);
fprintf('Distance between B-sat1 and Sat-2: %.2f kms\n', R23);
fprintf('Distance between Sat-2 and B-Sat2: %.2f kms\n', R34);