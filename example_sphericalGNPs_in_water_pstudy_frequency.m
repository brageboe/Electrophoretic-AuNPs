% 29/7/2024 Brage Bøe Svendsen
% Parameter study: background-subtracted heatrate wrt operating frequency
% Ambient: water w/ no conductivity
% Scatterer: GNP suspension
% Net charge: depends on GNP radius
clear variables
clc
%% PARAMETERS

% Constants
% eps0 = 8.8541878171e-12;   % vacuum permittivity
e0 = 1.602176634e-19;      % electron charge
% rho = 1.93e4;               % Au mass density

% System
E0 = 2.3e3;     % [V/m]
t = 300;        % exposure time [s]
T0 = 25;        % starting temperature [C]

% Sample
eta = 0.87e-3;      % viscosity of ambient/host
phi = 2e-2;         % volume fraction of GNPs
V_sample = 2e-18;    % Arbitrary, just be within electrostatic approx.
R_sample = (3*V_sample/(4*pi))^(1/3);
sigma_bg = 0;       % conductivity of ambient
R_GNP = 1.65e-9;    % radius of each GNP

%% FUNCTION HANDLES
q = @(R) net_charge_Rostalski(2*R,2*R);
% q = @(R) 6.7*e0; % used for testing, comparing to Collins sample B

eps_ambient = @(f) permittivity_water(f, sigma_bg);
eps_suspension = @(f, R) mean(permittivity_eph_suspension(f, phi, eta, eps_ambient(f), q(R), R*2, R*2)); % Spherical GNP, so take mean[eps_x eps_y eps_z] to get 1x1 size

C_abs = @(f, R) mean(C_abs_spheroid(f, eps_suspension(f, R), eps_ambient(f), R_sample*2, R_sample*2));
C_amb = @(f) C_amb_spheroid(f, eps_ambient(f), R_sample*2, R_sample*2);
F_abs = @(f, R) C_abs(f, R) / C_amb(f);

%% CALCULATE RESULTS
% Heatrate @ frequency sweep
freqs = logspace(6, 9, 1e3); 
T = zeros(1, length(freqs));
T_bg = zeros(1, length(freqs));
for i=1:length(freqs)
    T_bg(i) = max(bioheat_uniform_water(freqs(i), E0, t, V_sample, C_amb(freqs(i)), eps_ambient(freqs(i)), T0), [], 'all');
    T(i) = max(bioheat_uniform_water(freqs(i), E0, t, V_sample, C_abs(freqs(i), R_GNP), eps_ambient(freqs(i)), T0), [], 'all');
end
heatrate_sweep = (T-T_bg)/t;

%% FIGURE
figure()
lw = 1.5;
semilogx(freqs, heatrate_sweep, 'LineWidth', lw)
hold on
ylabel("$dT/dt$ [$^\circ$C/s]", 'interpreter', 'latex')
xlabel("$f$ [Hz]", 'interpreter', 'latex')
% title("Background-subtracted heat rate")
% subtitle("$(T-T_{bg})/t$",'interpreter','latex')
grid on
hold off

figure()
semilogx(freqs, T-T0, 'LineWidth', lw)
hold on
ylabel("$\Delta T$ [$^\circ$C]", 'interpreter', 'latex')
xlabel("$f$ [Hz]", 'interpreter', 'latex')
title("Temperature change after "+round(t/60)+" minutes")
subtitle("$T-T_{0}$",'interpreter','latex')
grid on
hold off

figure()
semilogx(freqs, T-T_bg, 'LineWidth', lw)
hold on
ylabel("$T$ [$^\circ$C]", 'interpreter', 'latex')
xlabel("$f$ [Hz]", 'interpreter', 'latex')
title("Background-subtracted temperature after "+round(t/60)+" minutes")
subtitle("$T-T_{bg}$",'interpreter','latex')
grid on
hold off

figure()
semilogx(freqs, T, 'LineWidth', lw)
hold on
semilogx(freqs, T_bg, 'LineWidth', lw)
ylabel("$T$ [$^\circ$C]", 'interpreter', 'latex')
xlabel("$f$ [Hz]",'interpreter','latex')
title("Max temperature after "+round(t/60)+" minutes")
legend("Suspension $T$","Ambient $T_{bg}$",'interpreter','latex')
grid on
hold off
%% FUNCTIONS
function T = bioheat_uniform_water(freq, E0, t, V_sample, C_abs, eps_host, T0) 
    % Return temperature T in suspension of volume V_sample
    c_w = 4186;     % Water specific heat capacity   [J/(kg*K)] 
    rho_w = 998;    % Water density                  [kg/m^3]   
    k_w = 0.6071;	% Water thermal conductivity     [W/(m*K)]  
    
    omega = 2*pi*freq;
    k0 = omega/299792458;
    mu0 = 1.25663706127e-6;     % vacuum permeability
    
    res = 3;
    T0_grid = T0.*ones(res, res, res);	% Initial temperature (3d matrix) [degC]
    S = ones(res, res, res);            % Heatrate due to external heat source (RF)

    k_host = k0*sqrt(eps_host);     
    I0 = (E0.*E0)*real(k_host)/(2*omega*mu0);
    
    w_abs = C_abs .* I0 / V_sample; 
    S = S .* w_abs/(rho_w*c_w); 
    
    region_side_length = V_sample ^ (1/3) * 0.01; 
    dx = region_side_length / res; 
    
    D = k_w/(rho_w*c_w);
    P = 0;  % blood perfusion rate
    Ta = 0; % arterial blood temperature

    T = bioheatExact(T0_grid, S, [D, P, Ta], dx, t);  
end
function Tsweep = bioheat_uniform_sweep(freq, E0space, t, V_sample, C_abs, C_amb, eps_host, T0)
    Tsweep = zeros(1,length(E0space));    
    for i=1:length(E0space)
        Ttmp = bioheat_uniform_water(freq, E0space(i), t, V_sample, C_abs, eps_host, T0);
        Ttmp = max(Ttmp(:));
        Ttmp_bg = bioheat_uniform_water(freq, E0space(i), t, V_sample, C_amb, eps_host, T0);
        Ttmp_bg = max(Ttmp_bg(:));
        Tsweep(i) = (Ttmp-Ttmp_bg)/t;
    end
end
function dTdt = heatrate_collins(freq, E0, eps_flow)
    c_w = 4186;     % Water specific heat capacity   [J/(kg*K)] 
    rho_w = 998;    % Water density                  [kg/m^3]   
    eps_0 = 8.8541878128e-12;
    omega = 2*pi*freq;
    
    dTdt = 0.5*omega*eps_0*imag(eps_flow)*E0^2/(rho_w*c_w);
end