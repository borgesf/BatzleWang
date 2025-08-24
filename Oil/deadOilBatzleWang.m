function [Rho_Oil, V_Oil, K_Oil] = deadOilBatzleWang(T, P, API)
% DEADOILBATZLEWANG Dead oil properties (density, velocity, bulk modulus) from Batzle–Wang.
%
% Implements the Batzle–Wang correlations for crude (dead) oil using:
%   • Eq. (14)  : density at surface from API
%   • Eq. (18)  : pressure correction to density
%   • Eq. (19)  : temperature correction to density
%   • Eq. (20b) : velocity as a function of API, T, and P
% Bulk modulus is computed from K = rho * V^2 (converted to GPa).
%
% Inputs:
%   T   - double array : Temperature (°C)
%   P   - double array : Pressure (MPa)
%   API - double array : API gravity
%
% Outputs:
%   Rho_Oil - double array : Oil density (kg/m^3)
%   V_Oil   - double array : Oil velocity (m/s)
%   K_Oil   - double array : Oil bulk modulus (GPa)
%
% Example:
%   [rho, v, k] = deadOilBatzleWang([60 70], [25 30], [35 35]);
%
% Written by Filipe Borges in August 2025

arguments (Input)
    T   (:,:) double {mustBeFinite}
    P   (:,:) double {mustBeFinite}
    API (:,:) double {mustBePositive}
end
arguments (Output)
    Rho_Oil (:,:) double
    V_Oil   (:,:) double
    K_Oil   (:,:) double
end

% -----------------------
% Sanity check on sizes
% -----------------------
if ~isequal(size(T), size(P), size(API))
    error('Inputs T, P, and API must be the same size.');
end

%% Density (start in g/cm^3, then convert to kg/m^3)

% Eq. (14): surface/reference density from API (rho0 in g/cm^3)
rho0 = 141.5 ./ (API + 131.5);

% Eq. (18): pressure-corrected density rho_P (still g/cm^3)
rhoP = rho0 ...
     + (0.00277 .* P - 1.71e-7 .* P.^3) .* (rho0 - 1.15).^2 ...
     + 3.49e-4 .* P;

% Eq. (19): temperature-corrected density (still g/cm^3)
rho_gcc = rhoP ./ (0.972 + 3.81e-4 .* (T + 17.78).^1.175);

% Convert to kg/m^3
Rho_Oil = 1000 .* rho_gcc;

%% Velocity (m/s)

% Eq. (20b): velocity in terms of API, T (°C) and P (MPa)
V_Oil = 15450 ./ sqrt(77.1 + API) ...
      - 3.77 .* T ...
      + 4.64 .* P ...
      + 0.0115 .* (0.36 .* sqrt(API) - 1) .* T .* P;

%% Bulk modulus (GPa)
% Using K = rho * V^2, with rho in kg/m^3 and V in m/s
K_Oil = Rho_Oil .* V_Oil.^2 * 1e-9;

end
