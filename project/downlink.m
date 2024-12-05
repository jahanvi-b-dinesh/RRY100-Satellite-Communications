clc;
close all;
c = 3e8; % Speed of light in meters per second
f_down = 13e9; % Downlink frequency in Hz
d = 42162.9e3; % Distance to GEO satellite in meters from the center of the Earth in meters
R = 39355.92e3; % Distance between ground station GOT and Satellite SAT1  
Rsquare = R^2;
G_max_tra = 56.19; %gain in dBi of ground station receiving antenna with diamater 8 meters and effeciency 0.7

%% 


% Satellite Antenna parameters
D = 3.2; % Diameter of the antenna in meters
eta = 0.7; % Efficiency of the antenna
P_t = 20.791812; %considering tranmistted power in 120 Watts converted in dBW (less than actual one)

%% 


%Downlink Performance characterstics
lamda_down = c / f_down;
G_max_down = (pi^2 * D^2 * eta) / lamda_down^2;
G_max_down_dB = 10 * log10(G_max_down); % Convert gain to decibels (dB)
fprintf('Downlink Antenna Gain: %f linear units\n', G_max_down);
fprintf('Downlink Antenna Gain: %f dBi\n', G_max_down_dB);

EIRP_max_down = P_t + G_max_down_dB;
fprintf('Downlink EIRP of Satellite: %f dBW\n', EIRP_max_down);

PFD_down = EIRP_max_down - 10*log10(4*pi*Rsquare); %Power Flux density dBW/m^2
fprintf('Downlink Power Flux Density: %f dBW/m^2\n', PFD_down);

L_FS_down = (4*pi*R / lamda_down)^2; %Attenuation of Free space
L_FS_down_dB = 10*log10(L_FS_down);
fprintf('Attenuation of free space: %f dBW/m^2\n', L_FS_down_dB);

P_R_down = EIRP_max_down - L_FS_down_dB + G_max_rec; %power received by Ground station antenna
fprintf('Power Received by ground station Antenna: %f dBW\n', P_R_down);

%% 

%Additional Losses
Fd_loss = 2; % Feeder loss typically ranges from 0.5 - 1 dB
Pol_loss = 1;  % Polarization mismatch loss is typically around 0.5 to 1 dB
Atm_loss = 1; % Atmospheric attenuation by atmospheric gases in dB are uaually 0.5 - 2 dB - assumption

theta_off = 0.1; % considering the offset angle as 0.1 degrees for calculating the pointing loss
theta_3dB_down = 70*lamda_down / D; % in degrees
A_point_loss = 12 * (theta_off / theta_3dB_down)^2; % Antenna Pointing loss in dB
fprintf('Antenna Pointing Loss: %f dB\n', A_point_loss);

% Parameters for rain attenuation calculation
R = 30; % Rain rate in mm/h as per (ITU-R P.837)
k = 0.03041; % k Coefficient for horizontal polarization at 13 GHz (ITU-R P.838-3)
alpha = 1.1586; % alpha Coefficient for horizontal polarization at 13 GHz (ITU-R P.838-3)

% Calculate specific attenuation - Rain Attenuation 
gamma_R = k * R^alpha;
% Display the result
fprintf('Specific Attenuation model for Rain (gamma_R): %f dB/km\n', gamma_R);

%Rain height model for prediction methods
h_o = 2.4; % isotherm height above mean sea level to be obtained from zip data of ITU-R P.839-4 

h_r = h_o + 0.36; % Rain height model as per ITU-R P.839-4 

Raintot_loss = h_r + gamma_R;
fprintf('total Rain loss: %f dBn\n', Raintot_loss);


%% 

P_R_down_tot = P_R_down - Fd_loss - Pol_loss - Atm_loss - A_point_loss - Raintot_loss; % Total power received by satellite antenna after considering all losses
fprintf('Power Received by Satellite Antenna after all losses: %f dBW\n', P_R_down_tot);

%% 

%Carrier Power to Noise Power spectral density (C/N_0)

L_ftx = 1.5; %feeder loss between amplifier and antenna in dB - assumption
L_frx = 1; %Loss between ground station antenna and reciever in dB - assumption
L_frx_con = 1.258; %converted from dB to real values for calculation of system noise temperature
T_F = 290; % Thermodynamic temperature of the connection in K - assumption based on pg 191
T_A = 65; % Antenna noise temperature in K - assumption 
T_eRX = 75; %Effective inpiut noise temperature which needs to be calculated by cascaded system in pg 191
L_r = 3; %page 216 reference in dB CHECK THIS ----------------------------- PAGE 191
Tm = 228.6; %Tranmission medium in dBW/HzK

T = T_A/L_frx_con + T_F*(1 - 1 / L_frx_con) + T_eRX; % system noise temperature 
fprintf('System Noise Temperature: %f K\n', T);

G_T = G_max_down_dB - L_r - L_frx - 10*log10(T);  % Figure of Mert of the satellite
fprintf('Figure of Merit (G/T) of Ground Station: %f dbK^-1\n', G_T);

C_No = EIRP_max_down - L_FS_down_dB + G_T + Tm;
fprintf('Downlink : Carrier Power to Noise power spectral density : %f dbHz\n', C_No);
%% Extended Code to Calculate Eb/No and BER

bandwidth_Hz = 13.8105e9;
beta = 0.3;
M = 4;
% Calculate SNR in dB
SNR_dB = C_No - 10 * log10(bandwidth_Hz);

% Display the result
fprintf('The SNR is %.2f dB\n', SNR_dB);

SpectralEfficiency = log2(M) / (1 + beta);

Eb_N0_dB = SNR_dB - 10 * log10(SpectralEfficiency);
fprintf('E_b/N_0: %.2f dB\n',Eb_N0_dB);

 


