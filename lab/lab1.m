clc;
close all;
clear all;

% Data from Table 1
Elevation = [90, 65, 55, 45, 40, 35, 30, 25]; % Elevation in degrees
m = [1, 1.103, 1.221, 1.414, 1.556, 1.743, 2, 2.366]; % Air mass factor
Vant = [0.8508, 0.8520, 0.8536, 0.8510, 0.8524, 0.8545, 0.8563, 0.8585]; % Antenna voltage (V)
Vwarm = [1.0891, 1.0954, 1.0944, 1.0919, 1.0931, 1.0934, 1.0947, 1.0945]; % Warm load voltage (V)
Vhot = [1.1319, 1.1375, 1.1362, 1.1349, 1.1356, 1.1362, 1.1373, 1.1370]; % Hot load voltage (V)
Twarm = [294, 294.1, 294.2, 294.3, 294.4, 294.5, 294.5, 294.6]; % Warm load temperature (K)
Thot = [345.9, 346.1, 346.4, 346, 346.3, 346.2, 346.2, 346]; % Hot load temperature (K)
Tground = [282.2, 282.2, 282.2, 282.2, 282.2, 282.2, 282.2, 282.2]; % Ground temperature (K)

% Constants
Teff = 0.95 * Tground;
Tbg = 2.8;
deltaThot = -1.8;

%calculate Trec
Y = Vhot / Vwarm;
Trec = (Thot - Y*Twarm)/(Y-1);
disp("Trec = ")
disp(Trec)

% Calculate antenna temperatures and m(epsilon)tauZ
V1 = Vwarm - Vant;
V2 = Vwarm - Vhot;
Tant = Twarm - (V1 ./ V2 .* (Twarm - (Thot + deltaThot)));
fprintf('Antenna Temperature: %.4f\n', Tant)
m_tauZ = -log((Teff - Tant) ./ (Teff - Tbg));

% Perform linear regression to find Tau_Z
p = polyfit(m, m_tauZ, 1);
tau_Z = p(1);  % Slope of the linear fit

% Display Tau_Z
fprintf('Final Tau_Z (tropospheric zenith opacity): %.4f\n', tau_Z);

% Plotting
figure;
plot(m, m_tauZ, '-o', 'DisplayName', 'm(epsilon)tauZ Data');
hold on;
x_fit = linspace(min(m), max(m), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '-', 'DisplayName', 'Linear Fit');
xlabel('m (Air mass factor)');
ylabel('m(epsilon)tau_Z');
title('Plot of m(epsilon)\tau_Z as a function of m(epsilon)');
grid on;
legend show;
line(xlim(), [0, 0], 'Color', 'red', 'LineStyle', '--'); % zero line for reference



% Calculate Tropospheric Transmission
transmission = exp(-tau_Z .* m);

% Calculate Tropospheric Attenuation
attenuation = 1 - transmission;

% Display the results
disp('Tropospheric Transmission:');
disp(transmission);

disp('Tropospheric Attenuation:');
disp(attenuation);