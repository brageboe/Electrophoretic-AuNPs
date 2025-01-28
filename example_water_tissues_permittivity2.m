% Permittivity of tissues and water.

% sigma_bg = [0,1e-3,1e-2, 5e-2, 1e-1, 1];       
% sigma_bg = [0,1e-3,3e-2,1e-1];    
sigma_bg = [0,1e-3,1e-1];    
freqs = logspace(5, 11, 1e3); 

lw = 2;
fs = 14;

% IMAGINARY PART
f1 = figure();
color = [0 0.4470 0.7410];
lc = length(sigma_bg)*5;
for i=1:length(sigma_bg)
    eps_ambient = @(f) permittivity_water(f, sigma_bg(i));
    loglog(freqs, imag(eps_ambient(freqs)), '-', 'LineWidth', lw, 'Color', color+[i/(lc/3) i/lc i/lc])
    hold on
end

eps_bfat = @(f) permittivity_tissue(2*pi*f, "Breast fat");
eps_btumor = @(f) permittivity_breastcancer(2*pi*f);
loglog(freqs, imag(eps_bfat(freqs)), '-', 'LineWidth', lw, 'Color', [0.9290 0.6940 0.1250])
loglog(freqs, imag(eps_btumor(freqs)), '-', 'LineWidth', lw, 'Color', [0.8500 0.3250 0.0980])
hold off

f1.Position = [680 458 560 300];
xlabel("Frequency (Hz)",'Interpreter','latex','FontSize',fs)
ylabel("$\varepsilon_r''$",'Interpreter','latex','FontSize',fs)
legend("$\sigma_{\mathrm{w}}=$ "+sigma_bg(1)+" S/m","$\sigma_{\mathrm{w}}=$ "+sigma_bg(2)+" S/m","$\sigma_{\mathrm{w}}=$ "+sigma_bg(3)+" S/m", ...
    "Breast fat", "Breast cancer", 'interpreter','latex','fontsize',fs-2, 'location', 'southeast')
ax1 = gca;
ax1.FontSize = fs+3;
ax1.Box = 'off';
ax1.Legend.NumColumns = 2;
ax1.TickLabelInterpreter = 'latex';
ax1.XLim = [min(freqs) max(freqs)];
ax1.XTick = [1e6 1e7 1e8 1e9 1e10];
ax1.XTickLabel = {"$10^6$" "" "$10^8$" "" "$10^{10}$"};
ax1.YLim = [1e-4 1e4];
ax1.YTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5];
ax1.YTickLabel = {"$10^{-4}$" "" "$10^{-2}$" "" "$10^{0}$" "" "$10^{2}$" "" "$10^{4}$"};
grid off

% REAL PART
f2 = figure();

eps_ambient = @(f) permittivity_water(f, 0);
loglog(freqs, real(eps_ambient(freqs)), '-', 'LineWidth', lw, 'Color', color)
hold on

eps_bfat = @(f) permittivity_tissue(2*pi*f, "Breast fat");
eps_btumor = @(f) permittivity_breastcancer(2*pi*f);
loglog(freqs, real(eps_bfat(freqs)), '-', 'LineWidth', lw, 'Color', [0.9290 0.6940 0.1250])
loglog(freqs, real(eps_btumor(freqs)), '-', 'LineWidth', lw, 'Color', [0.8500 0.3250 0.0980])
hold off

f2.Position = [680 458 560 300];
xlabel("Frequency (Hz)",'Interpreter','latex','FontSize',fs)
ylabel("$\varepsilon_r'$",'Interpreter','latex','FontSize',fs)
legend("Water", "Breast fat", "Breast cancer", 'interpreter','latex','fontsize',fs-2, 'location', 'southwest')
ax2 = gca;
ax2.FontSize = fs+3;
ax2.Box = 'off';
ax2.TickLabelInterpreter = 'latex';
ax2.XLim = [min(freqs) max(freqs)];
ax2.XTick = [1e6 1e7 1e8 1e9 1e10];
ax2.XTickLabel = {"$10^6$" "" "$10^8$" "" "$10^{10}$"};
ax2.YLim = [1e0 1e2];
grid off


