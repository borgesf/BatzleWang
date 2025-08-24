% BRINE_DEMO_SCRIPT
% Standalone script to demonstrate brineBatzleWang calculations.
%
% It generates three plots:
%   1) Velocity and density vs salinity (at fixed T and P)
%   2) Velocity and density vs temperature (at fixed P and S)
%   3) Velocity and density vs pressure (at fixed T and S)
%
% Written by Filipe Borges in August 2025

clear; clc; close all;

%% Fixed values
T_const = 60;      % °C
P_const = 30;      % MPa
S_const = 35000;   % ppm

%% 1) Velocity & density vs Salinity
S_vec = linspace(0, 180000, 150)';
T_S   = T_const * ones(size(S_vec));
P_S   = P_const * ones(size(S_vec));
[rho_S, v_S] = brineBatzleWang(T_S, P_S, S_vec);

figure('Position', [200, 200, 900, 520]);
yyaxis left
plot(S_vec, v_S, 'b-', 'LineWidth', 1.8); grid on; box on
ylabel('Velocity (m/s)')
yyaxis right
plot(S_vec, rho_S, 'r--', 'LineWidth', 1.8)
ylabel('Density (kg/m^3)')
xlabel('Salinity (ppm)')
title(sprintf('Velocity & Density vs Salinity (T = %.1f °C, P = %.1f MPa)', T_const, P_const))
legend('Velocity', 'Density', 'Location', 'best')

%% 2) Velocity & density vs Temperature
T_vec = linspace(0, 120, 150)';
P_T   = P_const * ones(size(T_vec));
S_T   = S_const * ones(size(T_vec));
[rho_T, v_T] = brineBatzleWang(T_vec, P_T, S_T);

figure('Position', [250, 250, 900, 520]);
yyaxis left
plot(T_vec, v_T, 'b-', 'LineWidth', 1.8); grid on; box on
ylabel('Velocity (m/s)')
yyaxis right
plot(T_vec, rho_T, 'r--', 'LineWidth', 1.8)
ylabel('Density (kg/m^3)')
xlabel('Temperature (°C)')
title(sprintf('Velocity & Density vs Temperature (P = %.1f MPa, S = %.0f ppm)', P_const, S_const))
legend('Velocity', 'Density', 'Location', 'best')

%% 3) Velocity & density vs Pressure
P_vec = linspace(0, 80, 150)';
T_P   = T_const * ones(size(P_vec));
S_P   = S_const * ones(size(P_vec));
[rho_P, v_P] = brineBatzleWang(T_P, P_vec, S_P);

figure('Position', [300, 300, 900, 520]);
yyaxis left
plot(P_vec, v_P, 'b-', 'LineWidth', 1.8); grid on; box on
ylabel('Velocity (m/s)')
yyaxis right
plot(P_vec, rho_P, 'r--', 'LineWidth', 1.8)
ylabel('Density (kg/m^3)')
xlabel('Pressure (MPa)')
title(sprintf('Velocity & Density vs Pressure (T = %.1f °C, S = %.0f ppm)', T_const, S_const))
legend('Velocity', 'Density', 'Location', 'best')
