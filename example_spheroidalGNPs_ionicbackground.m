clear variables
clc
%% %%%% PARAMETERS %%%% %%

% Constants
e0 = 1.602176634e-19;       % electron charge
resolution = 50; 

% System
freq = 10e6;
axis_interval = [1e-9 1e-6]; % oblate: D, prolate: L
tr_fixed = 3e-9;    % oblate: L, prolate: D
sigma_interval = [1e-4 1]; % sigma_w

E0 = 2.3e3;     % incident field strength [V/m]
t = 300;        % exposure time [s]ticks
T0 = 25;        % starting temperature [C]

% Sample
phi = 2e-2;     % volume fraction
R_susp = 5e-6;  % radius of suspension / cancer cell
tissue = "Water";

eta = 0.87e-3;
% sigma_bg = 0;           % ionic conductivity
eps_ambient = @(sigma_w) permittivity_water(freq, sigma_w);


%% %%%% FUNCTION HANDLES %%%% %%
q = @(D, L) net_charge_Rostalski(D, L);   % Sassaroli eq (25) using an effective spherical radius

eps_suspension = @(D, L, sigma_w) permittivity_eph_suspension(freq, phi, eta, eps_ambient(sigma_w), q(D,L), D, L, includeAmbientMD2017=true); % includeAmbientMD2017=true advised for tissues

C_abs = @(D, L, sigma_w) mean(C_abs_spheroid_orientation(freq, eps_suspension(D, L, sigma_w), eps_ambient(sigma_w), R_susp*2, R_susp*2, report_warning=true));
C_amb = @(sigma_w) C_amb_spheroid(freq, eps_ambient(sigma_w), R_susp*2, R_susp*2);
F_abs = @(D, L, sigma_w) C_abs(D, L, sigma_w) / C_amb(sigma_w);

T_bg = @(D, L, sigma_w) max(bioheat_uniform_water(freq, E0, t, R_susp, C_amb(sigma_w), eps_ambient(sigma_w), T0), [], 'all');
T = @(D, L, sigma_w) max(bioheat_uniform_water(freq, E0, t, R_susp, C_abs(D, L, sigma_w), eps_ambient(sigma_w), T0), [], 'all');

heatrate = @(D, L, sigma_w) (T(D, L, sigma_w) - T_bg(D, L, sigma_w))/t;
% heatrate = @(f, phi, D, L) T_bg(f, phi, D, L)/t;

friction_constant = @(eta, D,L) mean(stokes_drag(eta, D, L));

%% FIGURE F_abs %%
% PROLATE
f1 = figure();
make_plot_2d(@(sigma_w, L) F_abs(tr_fixed, L, sigma_w), sigma_interval, axis_interval, ...
    zlabel="$F_{abs}$", res=resolution);
% set_standard_plot_style();
% set_small_plot_style();
set_colorbar_plot_style();
% xlabel("$\sigma_{\textrm{w}}$ [S/m]");
ylabel("$L$ [m]");

title("Prolates ($D=$ "+tr_fixed*1e9+" nm)")
clim([1 538.8663]) 

ax1 = gca;
ax1.YTick = [1e-9 1e-8 1e-7 1e-6];
ax1.XTick = [1e-4 1e-3 1e-2 1e-1 1e0];
ax1.Position = [0.180484163826183,0.209466663734118,0.724515836173818,0.673566669432322];

% OBLATE
f2 = figure();
make_plot_2d(@(sigma_w, D) F_abs(D, tr_fixed, sigma_w), sigma_interval, axis_interval, ...
    zlabel="$F_{abs}$", res=resolution);
% set_standard_plot_style();
% set_small_plot_style();
set_colorbar_plot_style();
xlabel("$\sigma_{\textrm{w}}$ [S/m]");
ylabel("$D$ [m]");

title("Oblates ($L=$ "+tr_fixed*1e9+" nm)")
clim([1 538.8663]) 

ax2 = gca;
ax2.YTick = [1e-9 1e-8 1e-7 1e-6];
ax2.XTick = [1e-4 1e-3 1e-2 1e-1 1e0];
ax2.Position = [0.180484163826183,0.209466663734118,0.724515836173818,0.673566669432322];

% SPHERE
f3 = figure();
make_plot_2d(@(sigma_w, D) F_abs(D, D, sigma_w), sigma_interval, axis_interval, ...
    zlabel="$F_{abs}$", res=resolution);
% set_standard_plot_style();
% set_small_plot_style();
set_colorbar_plot_style();
% xlabel("$\sigma_{\textrm{w}}$ [S/m]");
ylabel("$D$ [m]");

title("Spheres ($D=L$)")
clim([1 538.8663]) 

ax3 = gca;
ax3.YTick = [1e-9 1e-8 1e-7 1e-6];
ax3.XTick = [1e-4 1e-3 1e-2 1e-1 1e0];
ax3.Position = [0.180484163826183,0.209466663734118,0.724515836173818,0.673566669432322];

%% FIGURE Temperature %%
% PROLATE
f4 = figure();
make_plot_2d(@(sigma_w, L) heatrate(tr_fixed, L, sigma_w), sigma_interval, axis_interval, ...
    zlabel="Heatrate $dT/dt$ [$^\circ$C/s]", res=resolution, isoline1=0.01, isoline2=0.001);
% set_standard_plot_style();
% set_small_plot_style();
set_colorbar_plot_style();
% xlabel("$\sigma_{\textrm{w}}$ [S/m]");
ylabel("$L$ [m]");

title("Prolates ($D=$ "+tr_fixed*1e9+" nm)")
clim([1e-4 0.049046])

ax4 = gca;
ax4.YTick = [1e-9 1e-8 1e-7 1e-6];
ax4.XTick = [1e-4 1e-3 1e-2 1e-1 1e0];
ax4.Position = [0.180484163826183,0.209466663734118,0.724515836173818,0.673566669432322];

% min = 0.00013316 @ [1, 1e-06]
% max = 0.036056 @ [0.0001, 1e-09]

% OBLATE
f5 = figure();
make_plot_2d(@(sigma_w, D) heatrate(D, tr_fixed, sigma_w), sigma_interval, axis_interval, ...
    zlabel="Heatrate $dT/dt$ [$^\circ$C/s]", res=resolution, isoline1=0.01, isoline2=0.001);
% set_standard_plot_style();
% set_small_plot_style();
set_colorbar_plot_style();
xlabel("$\sigma_{\textrm{w}}$ [S/m]");
ylabel("$D$ [m]");

title("Oblates ($L=$ "+tr_fixed*1e9+" nm)")
clim([1e-4 0.049046])

ax5 = gca;
ax5.YTick = [1e-9 1e-8 1e-7 1e-6];
ax5.XTick = [1e-4 1e-3 1e-2 1e-1 1e0];
ax5.Position = [0.180484163826183,0.209466663734118,0.724515836173818,0.673566669432322];

% min = 9.4865e-05 @ [1, 1e-06]
% max = 0.04106 @ [0.0001149, 1e-09]

% SPHERE
f6 = figure();
make_plot_2d(@(sigma_w, D) heatrate(D, D, sigma_w), sigma_interval, axis_interval, ...
    zlabel="Heatrate $dT/dt$ [$^\circ$C/s]", res=resolution, isoline1=0.01, isoline2=0.001);
% set_standard_plot_style();
% set_small_plot_style();
set_colorbar_plot_style();
% xlabel("$\sigma_{\textrm{w}}$ [S/m]");
ylabel("$D$ [m]");

title("Spheres ($D=L$)")
clim([1e-4 0.049046])

ax6 = gca;
ax6.YTick = [1e-9 1e-8 1e-7 1e-6];
ax6.XTick = [1e-4 1e-3 1e-2 1e-1 1e0];
ax6.Position = [0.180484163826183,0.209466663734118,0.724515836173818,0.673566669432322];

% min = 2.3108e-07 @ [1, 1e-09]
% max = 0.049046 @ [0.0001, 1e-06]