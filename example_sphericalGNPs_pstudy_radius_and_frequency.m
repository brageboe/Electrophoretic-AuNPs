%% 27/8/2024: not working
% plan was to merge two plots (dT/dt vs R, and dT/dt vs freq) into one
% figure with different x-axes. This is problematic when one axis is linear
% and the other logarithmic.

% clear variables
% clc
%% PARAMETERS

% Constants
% eps0 = 8.8541878171e-12;   % vacuum permittivity
e0 = 1.602176634e-19;      % electron charge
% rho = 1.93e4;               % Au mass density

% System
freq_const = 13.56e6;
E0 = 2.3e3;     % [V/m]
t = 300;        % exposure time [s]
T0 = 37;        % starting temperature [C]

% Sample
phi = 2e-2;
% V_sample = 2e-18;    % Arbitrary, just be within electrostatic approx.
% R_sample = (3*V_sample/(4*pi))^(1/3);
R_sample = 5e-6; V_sample = (4/3)*pi*R_sample^3;
sigma_bg = 0.0;%1e-3;       % conductivity of ambient
R_GNP_const = 1.65e-9;
tissue = "Breast fat";

% Setting ambient/host material properties
if tissue == "Breast fat"   
    % Viscosity of tissue, Sinkus2005 DOI: 10.1016/j.mri.2004.11.060
    eta = 0.55;     
    % Gabriel1996 http://niremf.ifac.cnr.it/docs/DIELECTRIC/AppendixC.html
    eps_ambient = @(f) permittivity_tissue(2*pi*f, "Breast fat"); 
elseif tissue == "Breast tumor"
    % Viscosity of tissue, Sinkus2005 DOI: 10.1016/j.mri.2004.11.060
    eta = 2.4;      
    % Breast cancer permittivity 0.5-20 GHz, Lazebnik DOI: 10.1088/0031-9155/52/20/002
    eps_ambient = @(f) permittivity_breastcancer(2*pi*f);   
elseif tissue == "Water"
    eta = 1e-3;
    sigma_bg = 0.03;           % ionic conductivity1e-7
    eps_ambient = @(f) permittivity_water(f, sigma_bg);
end

%% FUNCTION HANDLES
q = @(R) net_charge_Rostalski(2*R,2*R);
% q = @(R) 6.7*e0; % used for testing, comparing to Collins sample B

eps_suspension = @(f, R) mean(permittivity_eph_suspension(f, phi, eta, eps_ambient(f), q(R), R*2, R*2, includeAmbientMD2017=true)); % Spherical GNP, so take mean[eps_x eps_y eps_z] to get 1x1 size

% C_abs = @(f, R) C_abs_spheroid_orientation(f, eps_suspension(f, R), eps_ambient(f), R_sample*2, R_sample*2, field_oriented=true); 
% C_abs = @(f, R) C_abs_spheroid_orientation(f, eps_suspension(f, R), eps_ambient(f), R_sample*2, R_sample*2, flow_oriented=true, D_GNP=R*2, L_GNP=R*2, eta=eta); 
C_abs = @(f, R) mean(C_abs_spheroid(f, eps_suspension(f, R), eps_ambient(f), R_sample*2, R_sample*2, report_warning=true));
C_amb = @(f) C_amb_spheroid(f, eps_ambient(f), R_sample*2, R_sample*2);
F_abs = @(f, R) C_abs(f) / C_amb(f);

eps_suspension_constantcharge = @(f, R, charge) mean(permittivity_eph_suspension(f, phi, eta, eps_ambient(f), charge, R*2, R*2));
C_abs_constantcharge = @(f, R, charge) mean(C_abs_spheroid(f, eps_suspension_constantcharge(f, R, charge), eps_ambient(f), R_sample*2, R_sample*2));
F_abs_constantcharge = @(f, R, charge) C_abs_constantcharge(f, R, charge) / C_amb(f);

%% CALCULATE RESULTS
res = 1e3;          % resolution, #datapoints
% Heatrate @ R sweep
R_GNP = linspace(1.5e-9, 40e-9, res); 
T_R_sweep = zeros(1, length(R_GNP));
if tissue == "Water"
    T_bg = max(bioheat_uniform_water(freq_const, E0, t, R_sample, C_amb(freq_const), eps_ambient(freq_const), T0), [], 'all');
    for i=1:length(R_GNP)
        T_R_sweep(i) = max(bioheat_uniform_water(freq_const, E0, t, R_sample, C_abs(freq_const, R_GNP(i)), eps_ambient(freq_const), T0), [], 'all');
    end
elseif tissue == "Breast fat"
    T_bg = max(bioheat_uniform_breastfat(freq_const, E0, t, R_sample, C_amb(freq_const), eps_ambient(freq_const), T0), [], 'all');
    for i=1:length(R_GNP)
        T_R_sweep(i) = max(bioheat_uniform_breastfat(freq_const, E0, t, R_sample, C_abs(freq_const, R_GNP(i)), eps_ambient(freq_const), T0), [], 'all');
    end
elseif tissue == "Breast tumor"
    T_bg = max(bioheat_uniform_breasttumor(freq_const, E0, t, R_sample, C_amb(freq_const), eps_ambient(freq_const), T0), [], 'all');
    for i=1:length(R_GNP)
        T_R_sweep(i) = max(bioheat_uniform_breasttumor(freq_const, E0, t, R_sample, C_abs(freq_const, R_GNP(i)), eps_ambient(freq_const), T0), [], 'all');
    end
end
heatrate_R_sweep = (T_R_sweep-T_bg)/t;

% Heatrate @ frequency sweep
f_min = 6; f_max = 10;
freqs = logspace(f_min, f_max, res); 
T_f_sweep = zeros(1, length(freqs));
T_bg = zeros(1, length(freqs));
if tissue == "Water"
    for i=1:length(freqs)
        T_bg(i) = max(bioheat_uniform_water(freqs(i), E0, t, R_sample, C_amb(freqs(i)), eps_ambient(freqs(i)), T0), [], 'all');
        T_f_sweep(i) = max(bioheat_uniform_water(freqs(i), E0, t, R_sample, C_abs(freqs(i), R_GNP_const), eps_ambient(freqs(i)), T0), [], 'all');
    end
elseif tissue == "Breast fat"
    for i=1:length(freqs)
        T_bg(i) = max(bioheat_uniform_breastfat(freqs(i), E0, t, R_sample, C_amb(freqs(i)), eps_ambient(freqs(i)), T0), [], 'all');
        T_f_sweep(i) = max(bioheat_uniform_breastfat(freqs(i), E0, t, R_sample, C_abs(freqs(i), R_GNP_const), eps_ambient(freqs(i)), T0), [], 'all');
    end
elseif tissue == "Breast tumor"
    for i=1:length(freqs)
        T_bg(i) = max(bioheat_uniform_breasttumor(freqs(i), E0, t, R_sample, C_amb(freqs(i)), eps_ambient(freqs(i)), T0), [], 'all');
        T_f_sweep(i) = max(bioheat_uniform_breasttumor(freqs(i), E0, t, R_sample, C_abs(freqs(i), R_GNP_const), eps_ambient(freqs(i)), T0), [], 'all');
    end
end
heatrate_f_sweep = (T_f_sweep-T_bg)/t;

%% FIGURES
% % dT/dt vs R
% figure()
% lw = 1.5;
% hold on
% plot(R_GNP/1e-9, heatrate_R_sweep, 'LineWidth', lw)
% ylabel("$dT/dt$ [$^\circ$C/s]", 'interpreter', 'latex')
% xlabel("$R_{GNP}$ [nm]", 'interpreter', 'latex')
% title("Background-subtracted heat rate")
% grid on
% hold off
% 
% dT/dt vs f
% figure()
% lw = 1.5;
% plot(freqs, heatrate_f_sweep, 'LineWidth', lw)
% hold on
% ylabel("$dT/dt$ [$^\circ$C/s]", 'interpreter', 'latex')
% xlabel("$f$ [Hz]", 'interpreter', 'latex')
% title("Background-subtracted heat rate")
% subtitle("$(T-T_{bg})/t$",'interpreter','latex')
% grid on
% hold off

% % % % % % % %
lw = 1.5;
fs = 14;
% yt = [0.005 0.015 0.025];
yt = [0 0.01 0.02 0.03];
% yt = [0.005 0.01 0.015 0.02 0.025];


% dT/dt vs R 
f1 = figure();
s1 = plot(R_GNP*1e9*2, heatrate_R_sweep, 'LineWidth', lw);
hold on
ylabel("$dT/dt$ [$^\circ$C/s]", 'FontSize', fs, 'interpreter', 'latex')
ylabel("Heatrate [$^\circ$C/s]", 'FontSize', fs, 'interpreter', 'latex')
xlabel("AuNP diameter [nm]", 'FontSize', fs,  'interpreter', 'latex')
yticks(yt)
% yticklabels({'0.005', '', '0.015', '', '0.025'})
f1.Position = [680 458 560 280]; % [x0 y0 width height]
s1.LineWidth = lw;
% set(gca, 'FontSize', fs-2)
% set(gca, 'YGrid', 'on')
hold off

ax1 = gca;
set(ax1, 'FontSize', fs+5)
% xticks([0 20 40 60 80])
xticks(ax1, [0 20 40 60 80])

% dT/dt vs f
f2 = figure();
hold on
s2 = plot(freqs, heatrate_f_sweep);
% ylabel("$dT/dt$ [$^\circ$C/s]", 'FontSize', fs, 'interpreter', 'latex')
ylabel("Heatrate [$^\circ$C/s]", 'FontSize', fs, 'interpreter', 'latex')
xlabel("Frequency [Hz]", 'FontSize', fs, 'interpreter', 'latex')
% xlim([10^f_min 10^f_max])
% ylim([min(heatrate_f_sweep) max(heatrate_f_sweep)])
% yticks(yt)

f2.Position = f1.Position;
f2.CurrentAxes.XScale = 'log';
f2.CurrentAxes.XAxisLocation = 'bottom';
s2.Color = "#A2142F";
s2.LineWidth = lw+.5;
ax2 = gca;
set(ax2, 'FontSize', fs-2)
% % set(gca, 'YGrid', 'on')
% yticks(ax2, [0.01 0.02 0.03])
% yticks(ax2, yt)
% yticklabels(ax2, {"0", "0.01", "0.02", "0.03"})

% hold off
% f = gcf;
% f.Units = 'pixels';
% f.Position = ([100, 100, 1000, 1000]);

% semilogx(freqs, heatrate_f_sweep, 'LineWidth', lw)
% ax2 = axes('Position', ax1.Position, 'Color', "none");
% ax2.XAxisLocation = 'top';
% ax2.XLim = [min(freqs), max(freqs)];
% ax2.YTick = [];


% xlabel("$f$ [Hz]", 'interpreter', 'latex')

% Fig 9
if tissue=="Water"
    ylim([0 0.03])
    % set(gca,'YTickLabel',sprintf('%4.1e|',yt))
    ax2.YRuler.Exponent = -2;
end
set(ax2, 'FontSize', fs+5)
xticks(ax2, [1e6 1e7 1e8 1e9 1e10])
ax2.TickLabelInterpreter = 'latex';
