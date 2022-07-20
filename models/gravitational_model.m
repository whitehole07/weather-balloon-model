function g = gravitational_model(z, M_planet, R_planet)
%GRAVITATIONAL_MODEL Summary of this function goes here
%   Detailed explanation goes here
    g = (6.67430e-11 * M_planet) / ((z + R_planet)^2);
end

