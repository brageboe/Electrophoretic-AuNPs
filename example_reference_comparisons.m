clear variables
clc
%% DEFAULT PARAMETERS %%

% Constants
e0 = 1.602176634e-19;       % electron charge
eps0 =  8.8541878188e-12;

% System
freq = 10e6;
D = 5e-9;
L = D;

E0 = 2.3e3;     % incident field strength [V/m]
t = 300;        % exposure time [s]
T0 = 25;        % starting temperature [C]

% Sample
phi = 2e-2;     % volume fraction
R_susp = 5e-6;  % radius of suspension / cancer cell
ref = "chen2019";

%% REF DATA %%

if ref == "default"
    eta = 1e-3;
    sigma_bg = 0.0;   
    eps_ambient = permittivity_water(freq, sigma_bg);
elseif ref == "collins2018" % using approximate values
    D = 3.3e-9; L = D;
    phi = 2e-2;
    freq = 13.56e6;
    % ambient:
    eta = 0.87e-3;
    sigma_bg = 0;
    eps_ambient = permittivity_water(freq, sigma_bg);
elseif ref == "narasimh2020"
    D = 2e-9; L = D;
    phi = 2e-6;
    freq = 27e6;
    t = 600;
    T0 = 37;
    % PBS buffer:
    eta = 0.1;      % found at random website
    sigma_bg = 0.92; % @ 27 MHz [Narasimh2020->Gabriel1996]
    % eps_ambient = permittivity_water(freq, sigma_bg)
    eps_ambient = 111.5 + 1i*sigma_bg/(eps0*2*pi*freq);
    %
    E0 = field_strength_converter(200, eps_ambient, pi*120e-3*120e-3)
elseif ref == "bayford2022"
    D = 1.9e-9; L = D; % two samples: 1.9nm, 3.4 nm
    % phi = 1e-3; % unclear
    t = 30*60;
    freq = 2.6e9;
    % volume fraction
    V_GNP = (4*pi*D^3/3);
    m_GNP = V_GNP*19300;
    n = 0.06; % not sure. mg/mL = kg/m3
    N = n/m_GNP; % #GNPs
    phi = V_GNP * N;
    % PBS buffer:
    eta = 0.1;      % found at random website
    sigma_bg = 1; % @ "10% PBS" => around 1 S/m [DOI: 10.1088/1742-6596/687/1/012101, Fig4]
    eps_ambient = permittivity_water(freq, sigma_bg);
    E0 = field_strength_converter(5, eps_ambient, pi*2.9e-2*2.9e-2)
elseif ref == "chen2019"
    D=25e-9; L=50e-9;
    phi = 5e-5; 
    t = 10*60;
    freq = 13.56e9;
    eta = 1e-3; % DI water
    sigma_bg = 0; % DI water
    eps_ambient = permittivity_water(freq, sigma_bg);
    % instrument dimension unclear. assuming 10 x 10 cm area
    E0 = field_strength_converter(15, eps_ambient, pi*10e-2*10e-2) % 
else
    error("Input string not available.")
end


q = @(D, L) net_charge_Rostalski(D, L);   % Sassaroli eq (25) using an effective spherical radius
eps_suspension = @(f, phi, D, L) permittivity_eph_suspension(f, phi, eta, eps_ambient, q(D,L), D, L, includeAmbientMD2017=true); % includeAmbientMD2017=true advised for tissues

% C_abs = @(f, phi, D, L) C_abs_spheroid_orientation(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, field_oriented=true); 
% C_abs = @(f, phi, D, L) C_abs_spheroid_orientation(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, flow_oriented=true, D_GNP=D, L_GNP=L, eta=eta, report_warning=true); 
C_abs = @(f, phi, D, L) mean(C_abs_spheroid_orientation(f, eps_suspension(f, phi, D, L), eps_ambient, R_susp*2, R_susp*2, report_warning=true));
C_amb = @(f) C_amb_spheroid(f, eps_ambient, R_susp*2, R_susp*2);
F_abs = @(f, phi, D, L) C_abs(f, phi, D, L) / C_amb(f);


T_bg = @(f, phi) max(bioheat_uniform_water(f, E0, t, R_susp, C_amb(f), eps_ambient, T0), [], 'all');
T = @(f, phi, D, L) max(bioheat_uniform_water(f, E0, t, R_susp, C_abs(f, phi, D, L), eps_ambient, T0), [], 'all');

heatrate_bg = T_bg(freq, phi)/t;
heatrate_bgsub = (T(freq, phi, D, L) - T_bg(freq, phi))/t;
heatrate_total = T(freq, phi, D, L) / t;

fprintf("Heatrate, total (AuNPs+buffer) = "+heatrate_total+" degC/s \n")
fprintf("Heatrate, background (buffer) = "+heatrate_bg+" degC/s \n")    
fprintf("Heatrate, background-subtracted (AuNPs) = "+heatrate_bgsub+" degC/s \n")

%% Functions
function E0 = field_strength_converter(P, eps_m, area_suspension)
    imp0 = 376.730313412;       % free space impedance
    E0 = sqrt(2*imp0*P / ((area_suspension)*real(sqrt(eps_m))));
end
