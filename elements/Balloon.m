classdef Balloon < handle
    %BALOON Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        % Mass
        m_pay      % Payload mass            [kg]
        m_gas      % Gas mass                [kg]
        m_bal      % Balloon mass            [kg]

        m_dry      % Dry mass                [kg]
        m_tot      % Total mass              [kg]

        % Material
        t0         % Balloon thickness       [m]
        ten_str    % Tensile strength        [MPa]
        Vb         % Balloon (rigid) volume  [m^3]
        
        % Gas
        V0         % Initial gas volume      [m^3]
        R_star     % Specific gas constant   [J/(kgK)]

        % Drag
        Cd         % Drag coefficient        [-]

        % State
        equil_alt            % Equilibrium altitude    [m]
        burst_alt            % Burst altitude          [m]
        bursted = false      % If balloon bursted      [bool]

        exp_model  % Expansion model         [str]
    end
    
    methods
        function obj = Balloon(m_pay, rho_b, t0, Vb, ten_str, V0, R_star, Cd, exp_model)
            %BALOON Construct an instance of this class
            %   Detailed explanation goes here
            % Save needed quantities
            obj.m_pay = m_pay; obj.V0 = V0; obj.Cd = Cd;
            obj.R_star = R_star; obj.exp_model = exp_model;
            obj.t0 = t0; obj.Vb = Vb; obj.ten_str = ten_str;

            % Initial conditions
            [~, p0, T0] = atmospheric_model(0);  % Initial pressure and temperature  [hPa, °C]

            % Gas mass
            obj.m_gas = ((p0*1e2)*obj.V0) / (obj.R_star*(T0+273.15));  % Mass [kg]

            % Balloon mass
            r = ((3*Vb)/(4*pi))^(1/3);
            obj.m_bal = (4/3)*pi*rho_b * ((r+t0)^3 - r^3);

            % Dry mass
            obj.m_dry = obj.m_pay + obj.m_bal;
            
            % Total mass
            obj.m_tot = obj.m_gas + obj.m_dry;

            % Equilibrium altitude
            obj.equil_alt = obj.get_equilibrium_altitude();

            % Burst altitude
            obj.burst_alt = get_burst_altitude(obj);
        end
        
        function [V, rho, p] = gas_model(obj, z)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % Atmosphere temperature and pressure at h=z
            [~, patm, Tatm] = atmospheric_model(z); 
            
            switch obj.exp_model
                case "infinite"
                    % Approximation
                    p = patm;

                    % Gas density
                    rho = (p*1e2) / (obj.R_star * (Tatm+273.15));

                    % Gas volume
                    V = obj.m_gas / rho;
                case "maximum"
                    % Approximation
                    p = patm;

                    % Gas density
                    rho = (p*1e2) / (obj.R_star * (Tatm+273.15));

                    % Gas volume
                    V = obj.m_gas / rho;

                    if V >= obj.Vb
                        % Max volume reached
                        V = obj.Vb;
                        
                        % Gas density
                        rho = obj.m_gas / V;

                        % Gas pressure
                        p = rho * obj.R_star * (Tatm+273.15) * 1e-2;
                    end

                    % Tangential stress [MPa]
                    r = ((3*V)/(4*pi))^(1/3);
                    sigma_t = ((p-patm) * r)/(2 * obj.t0) * 1e-4;

                    if (sigma_t >= obj.ten_str) && ~obj.bursted
                        obj.bursted = true;
                    end
                case "constrained"
                    % Gas volume
                    V = obj.V0;

                    % Gas density
                    rho = obj.m_gas / V;

                    % Gas pressure
                    p = rho * obj.R_star * (Tatm+273.15) * 1e-2;
            end
        end
    
        function equil_alt = get_equilibrium_altitude(obj)
            function rho_gas = gas(z)
                [~, rho_gas, ~] = obj.gas_model(z);
            end
            function rho_air = air(z)
                [rho_air, ~, ] = atmospheric_model(z);
            end

            eq = @(z) air(z) - gas(z)*(obj.m_tot/obj.m_gas);
            equil_alt = fzero(eq, 0);
            
            % Reset state
            obj.bursted = 0;
        end

        function burst_alt = get_burst_altitude(obj)
            function p_gas = p_gas(z)
                [~, ~, p_gas] = obj.gas_model(z);
            end
            function r = r(z)
                [V_gas, ~, ~] = obj.gas_model(z);
                r = ((3*V_gas)/(4*pi))^(1/3);
            end
            function p_air = p_air(z)
                [~, p_air, ] = atmospheric_model(z);
            end

            eq = @(z) obj.ten_str - ((p_gas(z) - p_air(z)) * r(z))/(2 * obj.t0) * 1e-4;
            burst_alt = fzero(eq, 0);
            
            % Reset state
            obj.bursted = 0;
        end
    end
end
   

