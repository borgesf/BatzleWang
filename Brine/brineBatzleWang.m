function [Rho_Brine, V_Brine, K_Brine] = brineBatzleWang(T, P, Salinity)
% brineBatzleWang Computes brine properties (density, velocity, and bulk modulus) 
% for arrays of temperature, pressure, and salinity values.
%
% This function calculates brine properties using the Batzle-Wang empirical
% relationships. It accepts temperature (°C), pressure (MPa), and salinity
% (ppm) as arrays of matching size (1D or 2D), and returns the density
% (kg/m³), velocity (m/s), and bulk modulus (GPa) as arrays of the same size.
%
% Inputs:
%   T         - double array : Temperature in degrees Celsius
%   P         - double array : Pressure in MPa
%   Salinity  - double array : Salinity in parts per million (ppm)
%
% Outputs:
%   Rho_Brine - double array : Brine density in kg/m³
%   V_Brine   - double array : Brine velocity in m/s
%   K_Brine   - double array : Brine bulk modulus in GPa
%
% Example:
%   [rho, v, k] = brineBatzleWang([60 70], [25 30], [35000 36000]);
%
% Written by Filipe Borges in August 2025

arguments (Input)
    T (:,:) double {mustBeFinite}
    P (:,:) double {mustBeFinite}
    Salinity (:,:) double {mustBeNonnegative}
end

arguments (Output)
    Rho_Brine (:,:) double
    V_Brine   (:,:) double
    K_Brine   (:,:) double
end

% =======================
% Sanity Check for Sizes
% =======================
if ~isequal(size(T), size(P), size(Salinity))
    error('Inputs T, P, and Salinity must be the same size.');
end

% Convert salinity from ppm to mass fraction
S = Salinity / 1e6;

% Pure water density (Eq. 27a in the paper)
Rho_W = 1 + 1e-6 .* ( ...
        -80 .* T - 3.3 .* T.^2 + 0.00175 .* T.^3 + ...
         489 .* P - 2 .* T .* P + 0.016 .* T.^2 .* P - 1.3e-5 .* T.^3 .* P - ...
         0.333 .* P.^2 - 0.002 .* T .* P.^2 );

% Brine density (Eq. 27b in the paper)
Rho_Brine = Rho_W + S .* (0.668 + 0.44 .* S + 1e-6 .* ( ...
    300 .* P - 2400 .* P .* S + ...
    T .* (80 + 3 .* T - 3300 .* S - 13 .* P + 47 .* P .* S) ));

% Water velocity using matrix W (Table 1 in the paper)
W = [...
    1402.85      1.524         3.437e-3     -1.197e-5; 
    4.871       -0.0111        1.739e-4     -1.628e-6;
   -0.04783      2.747e-4     -2.135e-6      1.237e-8;
    1.487e-4    -6.503e-7     -1.455e-8      1.327e-10;
   -2.197e-7     7.987e-10     5.230e-11    -4.614e-13];

% Flatten arrays for element-wise computation
T_flat = T(:);
P_flat = P(:);
S_flat = S(:);

% Preallocate result vectors
V_Brine_flat = zeros(size(T_flat));

for i = 1:numel(T_flat)
    t = T_flat(i);
    p = P_flat(i);
    s = S_flat(i);

    % Calculate pure water velocity Vw (Eq. 28 in the paper)
    T_vec = [1, t, t^2, t^3, t^4];
    P_vec = [1; p; p^2; p^3];
    Vw = T_vec * (W * P_vec); 

    % Brine velocity at this element (Eq. 29 in the paper)
    V_Brine_flat(i) = Vw + ...
        s * (1170 - 9.6*t + 0.055*t^2 - 8.5e-5*t^3 + 2.6*p - 0.0029*t*p - 0.0476*p^2) + ...
        s^1.5 * (780 - 10*p + 0.16*p^2) - 820*s^2;
end

% Reshape to original input shape
V_Brine = reshape(V_Brine_flat, size(T));

% Converting density to kg/m³
Rho_Brine = 1000 * Rho_Brine;

% Bulk modulus calculation (in GPa)
K_Brine = (Rho_Brine) .* V_Brine.^2 * 1e-9;

end
