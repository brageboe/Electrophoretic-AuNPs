clear variables

tissue = "Breast cancer";

if tissue == "Breast cancer"
    eps_m = @(f, tissue) permittivity_breastcancer(2*pi*f);
else
    eps_m = @(f, tissue) permittivity_tissue(2*pi*f, tissue);
end

c_0 = 299792458;
k_imag = @(f, tissue) (2*pi*f/c_0) * sqrt(0.5*abs(eps_m(f, tissue) - 0.5*real(eps_m(f, tissue))));
k_real = @(f, tissue) (2*pi*f/c_0) * sqrt(0.5*abs(eps_m(f, tissue) + 0.5*real(eps_m(f, tissue))));

%% Attenuated part of E-field at depth r, normalized to amplitude E0
E = @(f, r, tissue) exp(-k_imag(f, tissue) * r );

%% Propagating part of E-field at depth r
% E = @(f, r, tissue) exp(1i * k_real(f, tissue) * r );

%% Combined, full, E-field at depth r
% E = @(f, r, tissue) real( exp(1i * k_real(f,tissue) * r) .* exp(-k_imag(f, tissue) * r) );

%% Result
% r = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5];
% f = 1000e6;
% E(f, r, "Breast fat")

r = 0.00:0.005:0.5;
fig = figure();
plot(r*100, E(1e6, r, tissue), 'LineWidth',1.5)
hold on
plot(r*100, E(10e6, r, tissue), 'LineWidth',1.5)
plot(r*100, E(100e6, r, tissue), 'LineWidth',1.5)
plot(r*100, E(1000e6, r, tissue), 'LineWidth',1.5)
ax = gca;
ax.FontSize = 14;
xlabel("Depth into tissue [cm]", 'interpreter','latex')
ylabel("$E/E_0$", 'interpreter','latex')
legend("1 MHz", "10 MHz", "100 MHz", "1 GHz", 'location', 'bestoutside')
xlim([min(r) max(r)].*100)
title("RF attenuation into "+tissue)
fig.Position = [300 300 600 300];

% ylim([-1 1])