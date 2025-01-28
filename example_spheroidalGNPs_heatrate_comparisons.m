% Same as 'example_surface_charge_comparisons.m' but without Gupta charge
clear variables 

%% %%%% PARAMETERS %%%% %%

% System
freq = 10e6;
A = linspace(3, 200, 1e3).*1e-9;  % length of one gnp axis (diameter or length)
B = min(A);     % length of other gnp axis (length or diameter)

E0 = 2.3e3;     % [V/m]
t = 300;        % exposure time [s]
T0 = 25;        % starting temperature [C]

% Sample
phi = 2e-2;     % volume fraction
R_susp = 5e-6;  % radius of suspension / cancer cell
tissue = "Water";

% ambient/host material properties
if tissue == "Breast fat"   
    % viscosity of tissue, Sinkus2005 DOI: 10.1016/j.mri.2004.11.060
    eta = 0.55;     
    % Gabriel1996 http://niremf.ifac.cnr.it/docs/DIELECTRIC/AppendixC.html
    eps_ambient = @(f) permittivity_tissue(2*pi*f, "Breast fat"); 
elseif tissue == "Breast tumor"
    % viscosity of tissue, Sinkus2005 DOI: 10.1016/j.mri.2004.11.060
    eta = 2.4;      
    % breast cancer permittivity 0.5-20 GHz, Lazebnik DOI: 10.1088/0031-9155/52/20/002
    eps_ambient = @(f) permittivity_breastcancer(2*pi*f);   
elseif tissue == "Water"
    eta = 0.87e-3;
    sigma_bg = 0.0;           % ionic conductivity
    eps_ambient = @(f) permittivity_water(f, sigma_bg);
else
    error("Input 'tissue' may only take value 'Breast fat', 'Breast tumor', or 'Water'.")
end

%% FUNCTION HANDLES
q = @(D, L) net_charge_Rostalski(D, L);   % Sassaroli eq (25) using an effective spherical radius

eps_suspension = @(f, phi, D, L) permittivity_eph_suspension(f, phi, eta, eps_ambient(f), q(D,L), D, L, includeAmbientMD2017=true); % [eps_transversal, eps_transversal, eps_axial]

% C_abs = @(f, phi, D, L) mean(C_abs_spheroid(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, approximation_check=true, report_warning=true));
C_amb = @(f) C_amb_spheroid(f, eps_ambient(f), R_susp*2, R_susp*2);

% C_abs = @(f, phi, D, L) C_abs_spheroid_orientation(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, field_oriented=true); 
C_abs = @(f, phi, D, L) C_abs_spheroid_orientation(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, flow_oriented=true, D_GNP=D, L_GNP=L, eta=eta, report_warning=true); 
% C_abs = @(f, phi, D, L) mean(C_abs_spheroid_orientation(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, report_warning=true, approximation_check=false));

if tissue == "Breast fat"   
    T_bg = @(f, phi, D, L) max(bioheat_uniform_breastfat(f, E0, t, R_susp, C_amb(f), eps_ambient(f), T0), [], 'all');
    T = @(f, phi, D, L) max(bioheat_uniform_breastfat(f, E0, t, R_susp, C_abs(f, phi, D, L), eps_ambient(f), T0), [], 'all');
elseif tissue == "Breast tumor"
    T_bg = @(f, phi, D, L) max(bioheat_uniform_breasttumor(f, E0, t, R_susp, C_amb(f), eps_ambient(f), T0), [], 'all');
    T = @(f, phi, D, L) max(bioheat_uniform_breasttumor(f, E0, t, R_susp, C_abs(f, phi, D, L), eps_ambient(f), T0), [], 'all');
elseif tissue == "Water"
    T_bg = @(f, phi, D, L) max(bioheat_uniform_water(f, E0, t, R_susp, C_amb(f), eps_ambient(f), T0), [], 'all');
    T = @(f, phi, D, L) max(bioheat_uniform_water(f, E0, t, R_susp, C_abs(f, phi, D, L), eps_ambient(f), T0), [], 'all');
end

heatrate = @(f, phi, D, L) (T(f, phi, D, L) - T_bg(f, phi, D, L))/t;
%%
heatrate_sphere = zeros(1,length(A));
heatrate_prolate = zeros(1,length(A));
heatrate_oblate = zeros(1,length(A));

for i=1:length(A)
    heatrate_sphere(i) = heatrate(freq, phi, A(i), A(i));
    heatrate_prolate(i) = heatrate(freq, phi, B, A(i));
    heatrate_oblate(i) = heatrate(freq, phi, A(i), B);
end
max_heatrate = max([max(heatrate_sphere), max(heatrate_prolate), max(heatrate_oblate)]);
min_heatrate = min([min(heatrate_sphere), min(heatrate_prolate), min(heatrate_oblate)]);

% SPHERES
f = figure();
subplot(1,3,1)
semilogy(A/1e-9, heatrate_sphere, '-')
hold on
title("Spheres")
xlabel("Diameter [nm]")
% ylabel("Heatrate $dT/dt$ [$^\circ$C/s]")
ylabel("Heatrate [$^\circ$C/s]")
ylim([min_heatrate max_heatrate])
grid on
hold off

% PROLATES
subplot(1,3,2)
semilogy(A/1e-9, heatrate_prolate, '--')
hold on
title("Prolates")
xlabel("Length [nm]")
ylim([min_heatrate max_heatrate])
grid on
hold off

% OBLATES
subplot(1,3,3)
semilogy(A/1e-9, heatrate_oblate, ':')
hold on
title("Oblates")
xlabel("Diameter [nm]")
ylim([min_heatrate max_heatrate])
grid on
hold off



axesHandles = get(gcf, 'children');
set(axesHandles, 'LineWidth', 1);
set(axesHandles, 'FontSize', 12);
set(axesHandles(1), 'YTickLabel', []);
set(axesHandles(2), 'YTickLabel', []);

f.Position = [680 458 560 280]; % [x0 y0 width height]
exportgraphics(f,"Figures/Fig7_heatrate_comparisons_minimum="+B*1e9+"nm.pdf")