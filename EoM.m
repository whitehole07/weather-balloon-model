function dxdt = EoM(~, x, balloon, M_planet, R_planet)
%EOM Summary of this function goes here
%   Detailed explanation goes here
% dxdt = [z_dot; z_dot_dot], x = [z; z_dot].

    if x(1) <= 0 && x(2) < 0
        dxdt = [0; 0];
        return
    end

    % State-dependant quantities
    [rho_air, ~, ~] = atmospheric_model(x(1));
    [V_gas, rho_gas, ~] = balloon.gas_model(x(1));
    g = gravitational_model(x(1), M_planet, R_planet);

    % Cross-sectional area
    r_sphere = ((3*V_gas)/(4*pi))^(1/3);
    Across = pi * r_sphere^2;

    % Compute contributions
    if ~balloon.bursted
        a = balloon.m_tot;
        b = 0.5 * rho_air * balloon.Cd * Across * sign(x(2));
        c = g * (balloon.m_dry + balloon.m_gas*(1 - (rho_air/rho_gas)));
    else
        a = balloon.m_dry;
        b = 0.5 * rho_air * balloon.Cd * Across * sign(x(2));
        c = g * balloon.m_dry;
    end

    % Compute output
    z_dot = x(2);
    z_dot_dot = -(b/a)*(x(2)^2) - c/a;
    
    % Assemble output
    dxdt = [z_dot; z_dot_dot];
end

