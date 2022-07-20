function [rho, p, T] = atmospheric_model(z)
%ISA Summary of this function goes here
%   Detailed explanation goes here
% Outputs:
% rho: air density      [kg/m^3]
% p  : air pressure     [hPa]
% T  : air temperature  [°C]

% Temperature and Pressure
if z <= 11000                      % Troposphere
    T = 15.04 - 0.00649*z;
    p = 101.29 * ((T+273.1)/288.08)^5.256;
elseif (z > 11000) && (z <= 25000) % Lower Stratosphere
    T = -56.46;
    p = 22.65 * exp(1.73 - 0.000157*z);
elseif z > 25000                   % Upper Stratosphere
    T = -131.21 + 0.00299*z;
    p = 2.488 * ((T+273.1)/216.6)^(-11.388);
end

% Density
rho = p / (0.2869 * (T + 273.1));

% Pression conversion (kPa -> hPa)
p = 10 * p;
end

