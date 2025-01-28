clear variables
clc
%% %%%% PARAMETERS %%%% %%

% Constants
e0 = 1.602176634e-19;       % electron charge

% System
freq = 10e6;
Dinterval = [1e-9 1e-6]; % Diameter
Linterval = [1e-9 1e-6]; % Length

E0 = 2.3e3;     % incident field strength [V/m]
t = 300;        % exposure time [s]
T0 = 25;        % starting temperature [C]

% Sample
phi = 2e-2;     % volume fraction
R_susp = 5e-6;  % radius of suspension / cancer cell
tissue = "Water";

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
    sigma_bg = 0.0;           % ionic conductivity
    eps_ambient = @(f) permittivity_water(f, sigma_bg);
else
    error("Input 'tissue' may only take value 'Breast fat', 'Breast tumor', or 'Water'.")
end

%% %%%% FUNCTION HANDLES %%%% %%
q = @(D, L) net_charge_Rostalski(D, L);   % Sassaroli eq (25) using an effective spherical radius
% q = @(D, L) net_charge_Gupta(D, L);         % Constant surface charge density

eps_suspension = @(f, phi, D, L) permittivity_eph_suspension(f, phi, eta, eps_ambient(f), q(D,L), D, L, includeAmbientMD2017=true); % includeAmbientMD2017=true advised for tissues

% C_abs = @(f, phi, D, L) C_abs_spheroid_orientation(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, field_oriented=true); 
% C_abs = @(f, phi, D, L) C_abs_spheroid_orientation(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, flow_oriented=true, D_GNP=D, L_GNP=L, eta=eta, report_warning=true); 
C_abs = @(f, phi, D, L) mean(C_abs_spheroid_orientation(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, report_warning=true));

% C_abs = @(f, phi, D, L) mean(C_abs_spheroid(f, eps_suspension(f, phi, D, L), eps_ambient(f), R_susp*2, R_susp*2, report_warning=true));
C_amb = @(f) C_amb_spheroid(f, eps_ambient(f), R_susp*2, R_susp*2);
F_abs = @(f, phi, D, L) C_abs(f, phi, D, L) / C_amb(f);

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
% heatrate = @(f, phi, D, L) T_bg(f, phi, D, L)/t;

% Needs further work before use: bioheat_grid
% V_susp = 4*pi*(R_susp^3)/3;
% T = @(f,phi,D,L) max(max(max(bioheat_grid(f, E0, t, V_susp, V_susp*10, eps_suspension(f, phi, D, L), eps_ambient(f), T0, 3, "Water"))), [], 'all');
% T_bg = @(f,phi,D,L) max(max(max(bioheat_grid(f, E0, t, V_susp, V_susp*10, eps_ambient(f), eps_ambient(f), T0, 3, "Water"))), [], 'all');
% heatrate = @(f, phi, D, L) (T(f, phi, D, L) - T_bg(f, phi, D, L))/t;
% heatrate = @(f, phi, D, L) T(f, phi, D, L)/t;

friction_constant = @(eta, D,L) mean(stokes_drag(eta, D, L));



%% %%%% FIGURES %%%% %%
% Uncomment sections to plot

%% FIGURE F_abs %%
% f = figure();
% make_plot_2d(@(D, L) F_abs(freq, phi, D, L), Dinterval, Linterval, ...
%     zlabel="$F_{abs}$", res=200, isoline1=1e3, isoline2=1);
% % set_standard_plot_style();
% set_small_plot_style();
% % set_colorbar_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% % title("$f=$ "+freq/1e6+" MHz; $\sigma_w=$ "+sigma_bg+" S/m")
% % title("$f=$ "+freq/1e6+" MHz; $\phi=$ "+phi+"; Ambient/host: "+tissue, "$\eta=$ "+eta+"; $\sigma_w=$"+sigma_bg+" S/m")
% % title("$\sigma_w=$ "+sigma_bg+" S/m")
% % % % % 
% % % settings for article Fig6 plots:
% title("$f=$ "+freq/1e6+" MHz; $\phi=$"+phi)
% clim([1 31045]) % both 
% ax1 = gca;
% ax1.YTick = [1e-9 1e-8 1e-7 1e-6];
% % % % %
% % % % settings for article Fig8 plots:
% % % ax1 = gca;
% % % clim([1 3102])    % max @ sigma_bg=0, eta=0.87e-3.
% % ax1.YTick = [1e-9 1e-8 1e-7 1e-6];

% f.Position = [9 2 6.500 4.3]; % short cbar plotstyle
% exportgraphics(f, "Figures/Fabs_DL_f="+freq/1e6+"MHz_phi="+phi+".pdf")

%% FIGURE Temperature %%
f2 = figure();
make_plot_2d(@(D,L) heatrate(freq, phi, D, L), Dinterval, Linterval, ...
    zlabel="Heatrate $dT/dt$ [$^\circ$C/s]", res=200, isoline1=0.01, minvalue=-1);
set_standard_plot_style();
% set_small_plot_style();
% set_colorbar_plot_style();
xlabel("$D$ [m]");
ylabel("$L$ [m]");
% title("Background-subtracted heat rate: $(T-T_{bg})/t$")
% title("$f=$ "+freq/1e6+" MHz; $\phi=$ "+phi+"; Ambient/host: "+tissue, "$\eta=$ "+eta);%+"; $\sigma_w=$"+sigma_bg+" S/m");
% title("$f=$ "+freq/1e6+" MHz; $\sigma_w=$ "+sigma_bg+" S/m")
% 
% % settings for article Fig7 plots:
title("$f=$ "+freq/1e6+" MHz")
% clim([2.0713e-07 0.041351])  % when eta=0.87e-3
clim([1e-4 0.041351])  % when eta=0.87e-3
hold on
min_r = 3e-9;
max_r = 100e-9;
color = [0 0.4470 0.7410 0.7];%"#0072BD"; % dark blue: [0 0.4470 0.7410]; light blue: [0.3010 0.7450 0.9330]. Fourth entry is alpha.
plot([min_r max_r], [min_r max_r], '-', 'Color', color); 
plot([min_r min_r], [min_r max_r], '--', 'Color', color); 
plot([min_r max_r], [min_r min_r], ':', 'Color', color);
ax1 = gca;
ax1.YTick = [1e-9 1e-8 1e-7 1e-6];
hold off
% exportgraphics(f2, "Figures/Fig7_heatrate_DL_f="+freq/1e6+"MHz_phi="+phi+".pdf")

% % settings for article Fig8 plots:
% title("$\sigma_w=$ "+sigma_bg+" S/m")
% ax2 = gca;
% clim([1e-6 0.041351])   % max @ sigma_bg=0, eta=0.87e-3
% % clim([1e-6 0.031262])   % max @ sigma_bg=1e-2
% ax2.YTick = [1e-9 1e-8 1e-7 1e-6];


%% FIGURE Charge %%
% f_charge = figure();
% make_plot_2d(@(D,L) q(D,L)/e0, Dinterval, Linterval, ...
%     zlabel="Net surface charge $q$ [$e_0$]", res=200, isoline1=10, isoline2=100, ticks=[1e1 1e2 1e3 1e4 1e5]);
% % set_standard_plot_style();
% set_small_plot_style();
% % set_colorbar_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% % f_charge.Position = [9 2 6.500 4.3];
% % exportgraphics(f_charge, "Figures/surfacecharge_DL_rostalski_small.pdf")

%% FIGURE Friction constant %%
% f_friction = figure();
% % make_plot_2d(@(D,L) friction_constant(eta, D, L), Dinterval, Linterval, zlabel="mean($\beta_i$) [m$\cdot$Pa$\cdot$s]", res=200, minvalue=1e-12)
% make_plot_2d(@(D,L) friction_constant(eta, D, L), Dinterval, Linterval, ...
%     zlabel="Average friction constant $\beta$ [mPa$\cdot$s]", res=200, ticks=[1e-11 1e-10 1e-9])
% % set_standard_plot_style();
% set_small_plot_style();
% % set_colorbar_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% % title("$\eta = $"+eta+" Pa$\cdot$s")
% % f_friction.Position = [1 3 6 5];
% % f_friction.Position = [9 2 6.500 4.3];
% % exportgraphics(f_friction, "Figures/friction_average_DL_eta="+eta+".pdf")


%% TESTING %%
% % Test: comparing k-wave with eq (11) in article2024, which % should give same result.
% figure()
% make_plot_2d(@(D,L) heatrate_nospatial_water_bgsub(freq, E0, R_susp, C_abs(freq,phi,D,L), C_amb(freq), eps_ambient(freq)), ...
%      Dinterval, Linterval, zlabel="$dT/dt$ (C$^\circ$/s)", res=200, isoline1=0.01)
% set_standard_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% title("$f = $"+freq/1e6+" MHz; $\phi = $"+phi+"; $\sigma_w=$"+sigma_bg+" S/m")
% hold off

% Test: Large background-subtracted heatrates even when Fabs=1.
% figure()
% make_plot_2d(@(D,L) C_abs(freq,phi,D,L), Dinterval, Linterval, zlabel="$C_{abs}", res=200, lim=-1)
% C_amb(freq)

function heatrate = heatrate_nospatial_water_bgsub(freq, E0, R_p, Cabs, Camb, host_permittivity)
    c_w = 4186;     % Water specific heat capacity   [J/(kg*K)] 
    rho_w = 998;    % Water density                  [kg/m^3]   
    omega = 2*pi*freq;
    k0 = omega/299792458;
    mu0 = 1.25663706127e-6;
    k_host = k0*sqrt(host_permittivity);     
    I0 = (E0.*E0)*real(k_host)/(2*omega*mu0);
    V_p = (4/3)*pi*R_p^3;

    heatrate = (Cabs - Camb) * I0 / (rho_w*c_w*V_p);
end



% % Test: does the condition \phi>\epsilon_{flow}/\epsilon_{host}$ lead to: 
% % -- F_abs < 1 (more absorption in ambient than suspension) ?
% % -- heatrate < 0 (negative background-subtracted heatrate) ?
% % under the assumptions:
% % -- F_abs and heatrate calculated with includeAmbient=true
% % (plot condition using includeAmbient=false to retrieve eps_flow)
% eps_flow = @(f, phi, D, L) permittivity_eph_suspension(f, phi, eta, eps_ambient(f), q(D,L), D, L, includeAmbient=false);
% % condition = @(D,L) real(mean(eps_flow(freq,phi,D,L))) / real(eps_ambient(freq));
% condition = @(D,L) imag(mean(eps_flow(freq,phi,D,L))) / imag(eps_ambient(freq));
% f_test = figure();
% make_plot_2d(@(D,L) condition(D,L), Dinterval, Linterval, zlabel="$\varepsilon_f/\varepsilon_h$", res=100, isoline1=phi, minvalue=1e-12);
% set_standard_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% title("Test: $\phi>\varepsilon_{flow}/\varepsilon_{host}$; $\phi=$ "+phi, "Ambient = "+tissue)
% 
% figure();
% make_plot_2d(@(D,L) real(mean(eps_suspension(freq,phi,D,L))), Dinterval, Linterval, zlabel="$\Re [\varepsilon_p]$", res=100, isoline1=phi, minvalue=1e-12);
% set_standard_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% title("$\varepsilon_p = (1-\phi)\varepsilon_m + \varepsilon_f$")
% 
% figure();
% make_plot_2d(@(D,L) imag(mean(eps_suspension(freq,phi,D,L))), Dinterval, Linterval, zlabel="$\Im [\varepsilon_p]$", res=100, isoline1=phi, minvalue=1e-12);
% set_standard_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% title("$\varepsilon_p = (1-\phi)\varepsilon_m + \varepsilon_f$")
% 
% eps_suspensionMD = @(f, phi, D, L) permittivity_eph_suspension(f, phi, eta, eps_ambient(f), q(D,L), D, L, includeAmbientMD2017=true); % includeAmbientMD2017=true advised for tissues
% figure();
% make_plot_2d(@(D,L) real(mean(eps_suspensionMD(freq,phi,D,L))), Dinterval, Linterval, zlabel="$\Re [\varepsilon_p]$", res=100, isoline1=phi, minvalue=1e-12);
% set_standard_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% title("$\varepsilon_p = \varepsilon_m + \varepsilon_f$")
% 
% figure();
% make_plot_2d(@(D,L) imag(mean(eps_suspensionMD(freq,phi,D,L))), Dinterval, Linterval, zlabel="$\Im [\varepsilon_p]$", res=100, isoline1=phi, minvalue=1e-12);
% set_standard_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% title("$\varepsilon_p = \varepsilon_m + \varepsilon_f$")
% 
% 
% ff1 = figure();
% make_plot_2d(@(D,L) abs(real(mean(eps_flow(freq,phi,D,L)))), Dinterval, Linterval, zlabel="$\Re [\varepsilon_f]$", res=100, isoline1=real(eps_ambient(freq)), minvalue=-1e2);
% set_standard_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% title("$\varepsilon_f$")
% 
% ff2 = figure();
% make_plot_2d(@(D,L) imag(mean(eps_flow(freq,phi,D,L))), Dinterval, Linterval, zlabel="$\Im [\varepsilon_f]$", res=100, isoline1=imag(eps_ambient(freq)), minvalue=1e-12);
% set_standard_plot_style();
% xlabel("$D$ [m]");
% ylabel("$L$ [m]");
% title("$\varepsilon_f$")