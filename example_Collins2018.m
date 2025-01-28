% 21/6/2024 Brage Bøe Svendsen
% Calculates heatrates based on the samples of Collins et al., 2018.
% (i)  Heatrate calculated as assumed done in Collins2018
% (ii) Heatrate calculated based on the suspension's absorption cross section.
% The goal is to 
% (1) Confirm that (i) results in the same values as in the plots in
% Collins2018 (Fig2 or FigS8, the latter is easier to read accurately)
% (2) Confirm that (ii) results are similar to (i)
%
% - Volume fractions are estimated from example_Collins2018_phi.m for any
% electric field strength E.
% - Buffer conductivities are set to zero: method (i) accounts only for
% imag(eps_flow) and ignores eps_ambient; method (ii) on the other hand 
% includes eps_ambient so we set its conductivity to zero to minimize the
% absorption from the ambient.
% - Volume of suspension V_sample can be arbitrary; it cancels out
% throughout the calculation. Ultimately, its only function is to check
% whether the domain is within the bounds of electrostatic approximation.
%
% Written by Brage Bøe Svendsen, 2024
clear variables
clc
%% PARAMETERS

% Constants
% eps0 = 8.8541878171e-12;   % vacuum permittivity
e0 = 1.602176634e-19;      % electron charge
% rho = 1.93e4;               % Au mass density

% System
freq = 13.56e6;
omega = 2*pi*freq;
E0 = 2.3e3;     % [V/m]
t = 300;        % exposure time [s]
T0 = 25;        % starting temperature [C]

% Sample
eta = 0.87e-3;      % buffer viscosity
R_GNP = 1.65e-9;
% See example_Collins2018_phi.m:
% phi = [0.0115, 0.0113, 0.0104]; % @ E0 = 3.0037e3
% phi = [0.0214, 0.0211, 0.0194]; % @ E0 = 2.2e3
phi = [0.0196, 0.0193, 0.0177]; % @ E0 = 2.3e3
V_sample = 2e-6;    % Sample volume, Collins2018: 2 mL. 
R_sample = (3*V_sample/(4*pi))^(1/3);
q = [-5.9*e0, -6.7*e0, -7.4*e0];        % net charge per GNP for samples [A,B,C]
sigma_bg = [4.1e-2, 1.2e-2, 1.16e-2];   % buffer conductivity for samples [A,B,C] unit [S/m]
% sigma_bg = [0,0,0];   

%% FUNCTION HANDLES
eps_ambient = @(f, sample) permittivity_water(f, sigma_bg(sample));
eps_suspension = @(f, sample) mean(permittivity_eph_suspension(freq, phi(sample), ...
    eta, eps_ambient(f, sample), q(sample), R_GNP*2, R_GNP*2, printContrast=false)); % Spherical GNP, so take mean[eps_x eps_y eps_z] to get 1x1 size

C_abs = @(f, sample) mean(C_abs_spheroid(f, eps_suspension(f, sample), eps_ambient(f, sample), R_sample*2, R_sample*2));
C_amb = @(f, sample) C_amb_spheroid(f, eps_ambient(f, sample), R_sample*2, R_sample*2);
F_abs = @(f, sample) C_abs(f, sample) / C_amb(f, sample);

%% CALCULATE RESULTS
% Heatrate @ E0 for all samples
for sample=1:3
    disp("%%%%%%% SAMPLE "+char(64+sample)+" %%%%%%%")
    eps_flow = mean(permittivity_eph_suspension(freq, phi(sample), eta, ...
        eps_ambient(freq, sample), q(sample), R_GNP*2, R_GNP*2, includeAmbient=false));
    heatrate_collins_bgsbtr = heatrate_collins(freq, E0, eps_flow)
    Fabs = F_abs(freq, sample)
    T_bg = bioheat_uniform(freq, E0, t, V_sample, C_amb(freq, sample), eps_ambient(freq, sample), T0);
    T_bg = max(T_bg(:))
    T = bioheat_uniform(freq, E0, t, V_sample, C_abs(freq, sample), eps_ambient(freq, sample), T0);
    T = max(T(:))
    heatrate = (T-T0)/t             % heating rate, including heating by background
    heatrate_bgsbtr = (T-T_bg)/t    % background subtracted heating rate
    relative_difference_bgsbtr = (heatrate_bgsbtr-heatrate_collins_bgsbtr)/heatrate_collins_bgsbtr
end

% Heatrate @ E0 sweep for all samples
E0sweep = linspace(1e2, 3e3, 500); 
Tsweep = zeros(3, length(E0sweep));
for sample=1:3
    Tsweep(sample,:) = bioheat_uniform_Esweep(freq, E0sweep, t, V_sample, ...
        C_abs(freq, sample), C_amb(freq, sample), eps_ambient(freq, sample), T0);
end
figure()
lw = 1.5;
hold on
plot(E0sweep, Tsweep, 'LineWidth', lw)
yline([0.0218, 0.0277, 0.0310], ':', 'LineWidth', lw-0.5)
xline(E0, ':', 'LineWidth', lw-0.5)
ylabel("$dT/dt$ [K/s]", 'interpreter', 'latex')
xlabel("$E$ [V/m]", 'interpreter', 'latex')
% title("Background-subtracted heatrate")
title({"$\phi=[$"+phi(1)+", "+phi(2)+", "+phi(3)+"$]$; $E_0=$"+E0+" V/m", ...
    "$\sigma_{water}=[$"+sigma_bg(1)+", "+sigma_bg(2)+", "+sigma_bg(3)+"$]$"}, 'interpreter', 'latex')
legend("A", "B", "C")
grid on
hold off

%% FUNCTIONS
function T = bioheat_uniform(freq, E0, t, V_sample, C_abs, eps_host, T0) 
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
%     T = w_abs/(rho_w*c_w);
end
function Tsweep = bioheat_uniform_Esweep(freq, E0space, t, V_sample, C_abs, C_amb, eps_host, T0)
    Tsweep = zeros(1,length(E0space));    
    for i=1:length(E0space)
        Ttmp = bioheat_uniform(freq, E0space(i), t, V_sample, C_abs, eps_host, T0);
        Ttmp = max(Ttmp(:));
        Ttmp_bg = bioheat_uniform(freq, E0space(i), t, V_sample, C_amb, eps_host, T0);
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