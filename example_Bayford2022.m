clear variables
% clc
%% Parameters
% constants
% eps0 = 8.8541878171e-12;   % vacuum permittivity
e0 = 1.602176634e-19;      % electron charge
rho = 1.93e4;               % Au mass density
% system
freq = 2.5999e9;
omega = 2*pi*freq;
E0 = 2e3; % [V/m]
t = 300;     % [s]
% sample
eta = 0.87e-3;
R_GNP = 1.65e-9;
%q = -5.9*e0;   % sample A
q = -6.7*e0;    % sample B
phi = 1e-2; 
V_sample = 2e-6; % Sample volume, Collins2018: 2 mL
R_sample = (3*V_sample/(4*pi))^(1/3);
T0 = 25;

%% Function handles
eps_ambient = @(f) permittivity_water(f, 1e-2);
eps_suspension = @(f) permittivity_eph_suspension(freq, phi, eta, eps_ambient(f), rho, q, R_GNP*2, R_GNP*2, printContrast=false);

C_abs = @(f) mean(C_abs_spheroid(f, eps_suspension(f), eps_ambient(f), R_sample*2, R_sample*2));
C_amb = @(f) C_amb_spheroid(f, eps_ambient(f), R_sample*2, R_sample*2);
F_abs = @(f) C_abs(f) / C_amb(f);

%% Calculate results
Fabs = F_abs(freq)
T = bioheat_uniform(freq, E0, t, V_sample, C_abs(freq), eps_ambient(freq), T0);
T = max(T(:))
heatrate = (T-T0)/t

%%
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
    q = ones(res, res, res);            % Volume rate of heat deposition (3d matrix) [W/m^3] of all ones

    k_host = k0*sqrt(eps_host);     
    I0 = E0*E0*real(k_host)/(2*omega*mu0);
    
    w_abs = C_abs .* I0 / V_sample; 
    q = q .* w_abs/(rho_w*c_w); 
    
    region_side_length = V_sample ^ (1/3) * 0.01; 
    dx = region_side_length / res; 
    
    D = k_w/(rho_w*c_w);
    P = 0;  % blood perfusion rate
    Ta = 0; % arterial blood temperature

    T = bioheatExact(T0_grid, q, [D, P, Ta], dx, t);  
end