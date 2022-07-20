%% Cleaning up
close all
clear
clc

%% Constants
% Planet
M_planet = 5.972e24;   % Planet gravitational constant [m^3/s^2]
R_planet = 6371010;    % Planet radius                 [m]

%% Variables
% Baloon
t0 = 10e-6;   % Balloon thickness [m]
Vb = 0.1131;  % Balloon maximum volume   [m^3]

% Gas
V0 = 0.11;         % Initial gas volume    [m^3]
R_star = 4124.2;   % Specific gas constant [J/(kgK)]

% Payload
m_pay = 0.110;   % Mass of the payload [kg]

% Material
rho_b = 1330;    % Balloon material density [kg/m^3]
ten_str = 34;    % Balloon tensile strength [MPa]

% Drag
Cd = 0.3;   % Drag coefficiennt [-]

% Expansion model
exp_model = "maximum";

%% Generate baloon
balloon = Balloon(m_pay, rho_b, t0, Vb, ten_str, V0, R_star, Cd, exp_model);

%% Solve equation of motion
% Initial state
z0 = 0;  % Initial altitude [m]
v0 = 0;  % Initial velocity [m/s]

% Condition
ToF = 0.5; % Time of flight [h]

%% Ode solver
% Event
opt = odeset('Events', @EoMStop);

% Solve
[t, x] = ode45(@(t, x) EoM(t, x, balloon, M_planet, R_planet), ...
               [0 ToF*3600], [z0; v0], opt);

%% Post-processing
% Unpack states
z = x(:, 1);    % Position [m]
v = x(:, 2);    % Velocity [m/s]

% Build rest if balloon fell
if t(end) < ToF*3600
    trest = (t(end)+1:1:ToF*3600)';
    t = [t; trest];
    z = [z; zeros(length(trest), 1)];
    v = [v; zeros(length(trest), 1)];
end

% Compute relevant quantities
W = []; B = []; D = []; rho_air = []; rho_He = []; V_He = []; dp = [];
sigma_t = []; 
for h = 1:length(z)
    gi = gravitational_model(z(h), M_planet, R_planet);
    [rho_airi, p_airi, ~] = atmospheric_model(z(h));
    [V_Hei, rho_Hei, p_Hei] = balloon.gas_model(z(h));

    rho_air = [rho_air; rho_airi];
    rho_He = [rho_He; rho_Hei];
    V_He = [V_He; V_Hei];
    dp = [dp; (p_Hei - p_airi)];

    r_sphere = ((3*V_Hei)/(4*pi))^(1/3);
    Across = pi * r_sphere^2;

    sigma_t = [sigma_t; ((p_Hei - p_airi) * r_sphere)/(2 * t0) * 1e-4];
    
    W = [W; balloon.m_tot*gi];
    B = [B; balloon.m_gas*gi*(rho_airi/rho_Hei)];
    D = [D; 0.5*rho_airi*balloon.Cd*Across];
end

%% Generate plots
% Altitude
figure("Position", [200 100 1500 900])
subplot(2, 3, 1)
plot(t./3600, z./1000, 'black', 'linewidth', 2)
grid on
hold on
if max(z./1000) >= 5; yline(12, 'r--', 'Tropopause', 'linewidth', 2); end
if max(z./1000) >= 40; yline(50, 'b--', 'Stratopause', 'linewidth', 2); end
if ~isnan(balloon.equil_alt)
yline(balloon.equil_alt/1000, '--', strcat(sprintf('%.1f', balloon.equil_alt/1000), ...
    'km - Equilibrium Altitude'), 'color', '#77AC30', 'linewidth', 2)
end
if ~isnan(balloon.burst_alt)
yline(balloon.burst_alt/1000, '--', strcat(sprintf('%.1f', balloon.burst_alt/1000), ...
    'km - Burst Altitude'), 'color', '#A2142F', 'linewidth', 2)
end
xlabel("Time [h]")
ylabel("Altitude [km]")
title("Altitude vs Time")

% Velocity
subplot(2, 3, 2)
plot(t./3600, v.*3.6, 'black', 'linewidth', 2)
grid on
hold on
yline(max(v.*3.6), '--', 'Max Velocity', 'color', '#D95319', 'linewidth', 2)
xlabel("Time [h]")
ylabel("Velocity [km/h]")
title("Velocity vs Time")

% Relative
subplot(2, 3, 3)
plot(z./1000, v.*3.6, 'black', 'linewidth', 2)
grid on
hold on
if max(z./1000) >= 5; xline(12, 'r--', 'Tropopause', 'linewidth', 2); end
if max(z./1000) >= 40; xline(50, 'b--', 'Stratopause', 'linewidth', 2); end
xlabel("Altitude [km]")
ylabel("Velocity [km/h]")
title("Velocity vs Altitude")

subplot(2, 3, 4)
plot(z./1000, W, 'color', '#A2142F', 'linewidth', 2)
hold on
plot(z./1000, B, 'color', '#77AC30', 'linewidth', 2)
plot(z./1000, D, 'color', '#EDB120', 'linewidth', 2)
grid on
legend("Weight", "Lift", "Drag")
hold on
xlabel("Altitude [km]")
ylabel("Force [N]")
title("Force vs Altitude")

subplot(2, 3, 5)
plot(z./1000, rho_air, 'color', '#4DBEEE', 'linewidth', 2)
hold on
plot(z./1000, rho_He, 'color', '#D95319', 'linewidth', 2)
legend("Air", "Helium")
grid on
hold on
xlabel("Altitude [km]")
ylabel("Density [kg/cu m]")
title("Density vs Altitude")

subplot(2, 3, 6)
plot(z./1000, V_He, 'color', '#7E2F8E', 'linewidth', 2)
grid on
hold on
xlabel("Altitude [km]")
ylabel("Volume [cu m]")
title("Baloon volume vs Altitude")

%% Balloon resistance
figure("Position", [200 100 1500 900])
subplot(2, 2, 1)
plot(t./3600, dp, 'color', '#EDB120', 'linewidth', 2)
grid on
hold on
xlabel("Time [h]")
ylabel("Pressure difference [hPa]")
title("Pressure difference vs Time")

subplot(2, 2, 2)
plot(t./3600, sigma_t, 'color', '#A2142F', 'linewidth', 2)
grid on
hold on
yline(ten_str, '--', 'Tensile Strength', 'color', '#7E2F8E', 'linewidth', 2)
xlabel("Time [h]")
ylabel("Tangential stress [MPa]")
title("Tangential stress vs Time")

subplot(2, 2, 3)
plot(z./1000, dp, 'color', '#EDB120', 'linewidth', 2)
grid on
hold on
xlabel("Altitude [km]")
ylabel("Pressure difference [hPa]")
title("Pressure difference vs Altitude")

subplot(2, 2, 4)
plot(z./1000, sigma_t, 'color', '#A2142F', 'linewidth', 2)
grid on
hold on
yline(ten_str, '--', 'Tensile Strength', 'color', '#7E2F8E', 'linewidth', 2)
xlabel("Altitude [km]")
ylabel("Tangential stress [MPa]")
title("Tangential stress vs Altitude")

