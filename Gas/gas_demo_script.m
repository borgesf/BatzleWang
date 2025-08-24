% GAS_DEMO_SCRIPT
% Standalone script to demonstrate gasBatzleWang calculations.
% It generates three plots where one variable varies and the other two are fixed:
%   1) Velocity and density vs Gas Gravity G
%   2) Velocity and density vs Temperature T
%   3) Velocity and density vs Pressure P
%
% Written by Filipe Borges in August 2025

clear; clc; close all;

%% Fixed (baseline) values
T_const = 60;     % degC
P_const = 25;     % MPa
G_const = 0.6;   % gas gravity (relative to air at 15.6°C, 1 atm)

%% 1) Velocity & density vs Gas Gravity (T, P fixed)
G_vec = linspace(0.55, 1.20, 150)';     % typical light to heavy gas range
T_G   = T_const * ones(size(G_vec));
P_G   = P_const * ones(size(G_vec));

[rho_G, v_G] = gasBatzleWang(T_G, P_G, G_vec);

figure('Position', [220, 220, 920, 520]);
yyaxis left
plot(G_vec, v_G, 'b-', 'LineWidth', 1.8); grid on; box on
ylabel('Velocity (m/s)')
yyaxis right
plot(G_vec, rho_G, 'r--', 'LineWidth', 1.8)
ylabel('Density (kg/m^3)')
xlabel('Gas gravity G (-)')
title(sprintf('Gas: Velocity & Density vs Gas Gravity (T = %.1f °C, P = %.1f MPa)', T_const, P_const))
legend('Velocity', 'Density', 'Location', 'best')

%% 2) Velocity & density vs Temperature (P, G fixed)
T_vec = linspace(0, 150, 150)';          % degC
P_T   = P_const * ones(size(T_vec));
G_T   = G_const * ones(size(T_vec));

[rho_T, v_T] = gasBatzleWang(T_vec, P_T, G_T);

figure('Position', [260, 260, 920, 520]);
yyaxis left
plot(T_vec, v_T, 'b-', 'LineWidth', 1.8); grid on; box on
ylabel('Velocity (m/s)')
yyaxis right
plot(T_vec, rho_T, 'r--', 'LineWidth', 1.8)
ylabel('Density (kg/m^3)')
xlabel('Temperature (°C)')
title(sprintf('Gas: Velocity & Density vs Temperature (P = %.1f MPa, G = %.2f)', P_const, G_const))
legend('Velocity', 'Density', 'Location', 'best')

%% 3) Velocity & density vs Pressure (T, G fixed)
P_vec = linspace(10, 80, 150)';         % MPa (start slightly above 0 to avoid divide-by-zero)
T_P   = T_const * ones(size(P_vec));
G_P   = G_const * ones(size(P_vec));

[rho_P, v_P] = gasBatzleWang(T_P, P_vec, G_P);

figure('Position', [300, 300, 920, 520]);
yyaxis left
plot(P_vec, v_P, 'b-', 'LineWidth', 1.8); grid on; box on
ylabel('Velocity (m/s)')
yyaxis right
plot(P_vec, rho_P, 'r--', 'LineWidth', 1.8)
ylabel('Density (kg/m^3)')
xlabel('Pressure (MPa)')
title(sprintf('Gas: Velocity & Density vs Pressure (T = %.1f °C, G = %.2f)', T_const, G_const))
legend('Velocity', 'Density', 'Location', 'best')
