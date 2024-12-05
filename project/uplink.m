clc;
close all;
c = 3e8; % Speed of light in meters per second
f_up = 14.4e9; % Uplink frequency in Hz 
d = 42162.9e3; % Distance to GEO satellite in meters from the center of the Earth in meters
R = 39355.92e3; % Distance between ground station GOT and Satellite BSAT2 ---- Wosrt case  
Rsquare = R^2;
G_max_rec = 51.87; %gain in dBi of Satellite receiving antenna with diamater 3.2 meters and effeciency 0.7

%% 

% Antenna parameters - Ground Station
D = 8; % Diameter of the antenna in meters
eta = 0.7; % Efficiency of the antenna
P_t = 21.76; %considering tranmistted power in 150 Watts converted in dBW
%P_t = 23.01; % 200W
%% 


%Uplink Performance characterstics
lamda_up = c / f_up;
G_max_up = (pi^2 * D^2 * eta) / lamda_up^2;
G_max_up_dB = 10 * log10(G_max_up); % Convert gain to decibels (dB)
fprintf('Uplink Antenna Gain: %f linear units\n', G_max_up);
fprintf('Uplink Antenna Gain: %f dBi\n', G_max_up_dB);

EIRP_max_up = P_t + G_max_up_dB;
fprintf('Uplink EIRP of Earth Station: %f dBW\n', EIRP_max_up);

PFD_up = EIRP_max_up - 10*log10(4*pi*Rsquare); %Power Flux density dBW/m^2
fprintf('Uplink Power Flux Density: %f dBW/m^2\n', PFD_up);

L_FS_up = (4*pi*R / lamda_up)^2; %Attenuation of Free space
L_FS_up_dB = 10*log10(L_FS_up);
fprintf('Attenuation of free space: %f dBW/m^2\n', L_FS_up_dB);

P_R_up = EIRP_max_up - L_FS_up_dB + G_max_rec; %power received by satellite antenna
fprintf('Power Received by Satellite Antenna: %f dBW\n', P_R_up);

%% 

%Additional Losses
Fd_loss = 2; % Feeder loss typically ranges from 1 to 3 dB - assumption

Pol_loss = 1;  % Polarization mismatch loss is typically around 0.5 to 1 dB - assumption

Atm_loss = 1; % Atmospheric attenuation by atmospheric gases in dB is usually around 0.5 - 2 dB - assumption

theta_off = 0.1; % considering the offset angle as 0.1 degrees for calculating the pointing loss
theta_3dB_up = 70*lamda_up / D; % in degrees
A_point_loss = 12 * (theta_off / theta_3dB_up)^2; % Antenna Pointing loss in dB -- Page 175 in m ytextbook
fprintf('Antenna Pointing Loss: %f dB\n', A_point_loss);

% Parameters for rain attenuation calculation
R = 30; % Rain rate in mm/h as per (ITU-R P.837)
k = 0.03738; % k Coefficient for horizontal polarization at 14 GHz (ITU-R P.838-3)
alpha = 1.1396; % alpha Coefficient for horizontal polarization at 14 GHz (ITU-R P.838-3)

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

P_R_up_tot = P_R_up - Fd_loss - Pol_loss - Atm_loss - A_point_loss - Raintot_loss ; % Total power received by satellite antenna after considering all losses
fprintf('Power Received by Satellite Antenna after all losses: %f dBW\n', P_R_up_tot);

%fprintf("Other losses: ", - Pol_loss - Atm_loss - A_point_loss - Raintot_loss);

%% 

%Carrier Power to Noise Power spectral density (C/N_0)

L_ftx = 1.5; %feeder loss between amplifier and antenna in dB - assumption
L_frx = 2; %Loss between satellite antenna and reciever in dB - assumption
L_frx_con = 1.58; %converted from dB to real values for calculation of system noise temperature
T_F = 290; % Thermodynamic temperature of the connection in K - assumption ----- Page 164
T_A = 290; % Antenna noise temperature in K - assumption --- Page 164
T_eRX = 290; %Effective inpiut noise temperature in K which needs to be calculated by cascaded system in pg 178
L_r = 3; %page 189 reference in dB
Tm = 228.6; %Tranmission medium in dBW/HzK (page 189 in my book)

T = T_A/L_frx_con + T_F*(1 - 1 / L_frx_con) + T_eRX; % system noise temperature (page 185 in my book)
fprintf('System Noise Temperature: %f K\n', T);

G_T = G_max_rec - L_r - L_frx - Pol_loss - 10*log10(T);  % Figure of Mert of the satellite (page 189 equation 5.40)
fprintf('Figure of Merit (G/T) of Satellite: %f dbK^-1\n', G_T);

C_No = EIRP_max_up - L_FS_up_dB + G_T + Tm; %page 188-189
fprintf('Uplink : Carrier Power to Noise power spectral density : %f dbHz\n', C_No);

%% Extended Code to Calculate Eb/No and BER

% % Required bandwidth and bit rate calculation
% bandwidth = 37.26e9; % Bandwidth in Hz (41.31 GHz)
% spectral_efficiency = 2; % QPSK has 2 bits per symbol
% R_b = bandwidth * spectral_efficiency; % Bit rate in bits per second
% fprintf("Bit Rate: %f b/s ", R_b)
% 
% % Calculate Eb/No from C/No
% Eb_No = C_No - 10*log10(R_b); % Eb/No in dB
% fprintf('Energy per Bit to Noise Power Density Ratio (Eb/No): %f dB\n', Eb_No);
% 
% % Convert Eb/No from dB to linear scale
% Eb_No_linear = 10^(Eb_No / 10);
% 
% % Theoretical BER for QPSK in AWGN
% ber_qpsk = 0.5 * erfc(sqrt(Eb_No_linear)); % BER calculation
% fprintf('Theoretical BER for QPSK: %e\n',ber_qpsk);


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



