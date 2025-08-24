% DEAD_OIL_DEMO_SCRIPT
% Standalone script to demonstrate deadOilBatzleWang calculations.
% It generates three plots where one variable varies and the other two are fixed:
%   1) Velocity and density vs API
%   2) Velocity and density vs Temperature
%   3) Velocity and density vs Pressure
%
% Written by Filipe Borges in August 2025

clear; clc; close all;

%% Fixed (baseline) values
T_const  = 60;     % degC
P_const  = 30;     % MPa
API_const = 35;    % API gravity

%% 1) Velocity & density vs API (T, P fixed)
API_vec = linspace(10, 50, 150)';          % API gravity range
T_API   = T_const * ones(size(API_vec));
P_API   = P_const * ones(size(API_vec));

[rho_API, v_API] = deadOilBatzleWang(T_API, P_API, API_vec);

figure('Position', [200, 200, 900, 520]);
yyaxis left
plot(API_vec, v_API, 'b-', 'LineWidth', 1.8); grid on; box on
ylabel('Velocity (m/s)')
yyaxis right
plot(API_vec, rho_API, 'r--', 'LineWidth', 1.8)
ylabel('Density (kg/m^3)')
xlabel('API gravity (°API)')
title(sprintf('Dead Oil: Velocity & Density vs API (T = %.1f °C, P = %.1f MPa)', T_const, P_const))
legend('Velocity', 'Density', 'Location', 'best')

%% 2) Velocity & density vs Temperature (P, API fixed)
T_vec = linspace(0, 120, 150)';             % degC
P_T   = P_const * ones(size(T_vec));
API_T = API_const * ones(size(T_vec));

[rho_T, v_T] = deadOilBatzleWang(T_vec, P_T, API_T);

figure('Position', [250, 250, 900, 520]);
yyaxis left
plot(T_vec, v_T, 'b-', 'LineWidth', 1.8); grid on; box on
ylabel('Velocity (m/s)')
yyaxis right
plot(T_vec, rho_T, 'r--', 'LineWidth', 1.8)
ylabel('Density (kg/m^3)')
xlabel('Temperature (°C)')
title(sprintf('Dead Oil: Velocity & Density vs Temperature (P = %.1f MPa, API = %.1f)', P_const, API_const))
legend('Velocity', 'Density', 'Location', 'best')

%% 3) Velocity & density vs Pressure (T, API fixed)
P_vec = linspace(0, 80, 150)';              % MPa
T_P   = T_const * ones(size(P_vec));
API_P = API_const * ones(size(P_vec));

[rho_P, v_P] = deadOilBatzleWang(T_P, P_vec, API_P);

figure('Position', [300, 300, 900, 520]);
yyaxis left
plot(P_vec, v_P, 'b-', 'LineWidth', 1.8); grid on; box on
ylabel('Velocity (m/s)')
yyaxis right
plot(P_vec, rho_P, 'r--', 'LineWidth', 1.8)
ylabel('Density (kg/m^3)')
xlabel('Pressure (MPa)')
title(sprintf('Dead Oil: Velocity & Density vs Pressure (T = %.1f °C, API = %.1f)', T_const, API_const))
legend('Velocity', 'Density', 'Location', 'best')
