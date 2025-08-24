function [Rho_Gas, V_Gas, K_Gas] = gasBatzleWang(T, P, G)
% GASBATZLEWANG Gas density, velocity, and bulk modulus from Batzle–Wang.
%
% Implements the Batzle–Wang correlations for hydrocarbon gases using
% pseudoreduced variables, compressibility factor Z, and an adiabatic
% correction for the bulk modulus.
%
% Referenced equations 
%   • Eq. (9a)–(9b): Pseudoreduced pressure P_pr and temperature T_pr
%   • Eq. (10a):     Real-gas density using compressibility Z
%   • Eq. (10b)–(10c): Z(P_pr, T_pr) and exponential term E
%   • Eq. (11b):     Approximation for adiabatic index γ0
%   • Eq. (11a):     Adiabatic bulk modulus K_s
%   • Velocity:      V = sqrt(K_s / ρ)
%
% Inputs:
%   T - double array : Temperature (°C)
%   P - double array : Pressure (MPa)
%   G - double array : Gas gravity (specific gravity to air at 15.6°C, 1 atm)
%
% Outputs:
%   Rho_Gas - double array : Gas density (kg/m^3)
%   V_Gas   - double array : Acoustic velocity (m/s)
%   K_Gas   - double array : Adiabatic bulk modulus (GPa)
%
% Example:
%   [rho, v, K] = gasBatzleWang([50 80], [20 40], [0.65 0.8]);
%
% Written by Filipe Borges in August 2025

arguments (Input)
    T (:,:) double {mustBeFinite}  % °C
    P (:,:) double {mustBeFinite}  % MPa
    G (:,:) double {mustBePositive} % gas gravity (relative to air)
end
arguments (Output)
    Rho_Gas (:,:) double
    V_Gas   (:,:) double
    K_Gas   (:,:) double
end

% -----------------------
% Sanity check on sizes
% -----------------------
if ~isequal(size(T), size(P), size(G))
    error('Inputs T, P, and G must have the same size.');
end

% Constants and helpful variables
T_a = T + 273.15;                        % absolute temperature (K)
R_MPa = 8.314e-3;                        % universal gas constant in MPa·m^3/(kmol·K)
M = 28.8 .* G;                           % molar mass (kg/kmol), air ≈ 28.8 → Eq. (2) logic

% ----------------------------------------------------------
% Eq. (9a)–(9b): Pseudoreduced variables (P_pr, T_pr)
% Pc ≈ 4.892 − 0.4048 G   [MPa]
% Tc ≈ 94.72 + 170.75 G   [K]
% ----------------------------------------------------------
Pc = 4.892 - 0.4048 .* G;                % MPa
Tc = 94.72 + 170.75 .* G;                % K
P_pr = P ./ Pc;
T_pr = T_a ./ Tc;

% ----------------------------------------------------------
% Eq. (10c): E(P_pr, T_pr)
% E = 0.109 (3.85 − T_pr)^2 * exp( −[0.45 + 8(0.56 − 1/T_pr)^2] * P_pr^1.2 / T_pr )
% ----------------------------------------------------------
exp_factor = (0.45 + 8 .* (0.56 - 1 ./ T_pr).^2) .* (P_pr .^ 1.2) ./ T_pr;
E = 0.109 .* (3.85 - T_pr).^2 .* exp(-exp_factor);

% ----------------------------------------------------------
% Eq. (10b): Z linear part in P_pr plus temperature terms
% Z = [0.03 + 0.00527 (3.5 − T_pr)^3] P_pr + (0.642 T_pr − 0.007 T_pr^4 − 0.52) + E
% ----------------------------------------------------------
Z = (0.03 + 0.00527 .* (3.5 - T_pr).^3) .* P_pr + ...
    (0.642 .* T_pr - 0.007 .* T_pr.^4 - 0.52) + E;

% Also compute ∂Z/∂P_pr at constant T_pr (needed in Eq. 11a)
% dZ/dP_pr = [0.03 + 0.00527 (3.5 − T_pr)^3] + dE/dP_pr
dE_dPpr = E .* ( - (0.45 + 8 .* (0.56 - 1 ./ T_pr).^2) .* (1.2 .* P_pr.^0.2) ./ T_pr );
dZ_dPpr = (0.03 + 0.00527 .* (3.5 - T_pr).^3) + dE_dPpr;

% ----------------------------------------------------------
% Eq. (10a): Real-gas density
% ρ = (M P) / (Z R T_a)
% Units here: M [kg/kmol], P [MPa], R [MPa·m^3/(kmol·K)], T_a [K]
% → ρ in kg/m^3
% ----------------------------------------------------------
Rho_Gas = (M .* P) ./ (Z .* R_MPa .* T_a);   % kg/m^3

% ----------------------------------------------------------
% Eq. (11b): Approximation for adiabatic index γ0 (dimensionless)
% γ0 ≈ 0.85 + 0.56/(P_pr + 2) + 27.1/(T_pr + 3.5)^2 − 8 * exp[ −0.65 (P_pr + 1) ]
% (Form consistent with the printed equation.)
% ----------------------------------------------------------
gamma0 = 0.85 + 0.56 ./ (P_pr + 2) + 27.1 ./ (T_pr + 3.5).^2 ...
         - 8 .* exp(-0.65 .* (P_pr + 1));

% ----------------------------------------------------------
% Eq. (11a): Adiabatic bulk modulus
% K_s = γ0 * P / [ 1 − (P_pr / Z) * (∂Z/∂P_pr)_T ]
% K_s in MPa → convert to GPa
% ----------------------------------------------------------
Ks_MPa = (gamma0 .* P) ./ ( 1 - (P_pr ./ Z) .* dZ_dPpr );
K_Gas  = Ks_MPa * 1e-3;                   % GPa

% Velocity from modulus and density
V_Gas = sqrt( (K_Gas * 1e9) ./ Rho_Gas ); % m/s   (convert GPa → Pa)

end
