% Permittivity of tissues and water.

% sigma_bg = [0,1e-3,1e-2, 5e-2, 1e-1, 1];       
sigma_bg = [0,1e-3,3e-2, 1e-1];       
freqs = logspace(0, 10, 1e3); 

figure();
lw = 1.5;
fs = 12;
color = [0 0.4470 0.7410];
lc = length(sigma_bg)*5;
for i=1:length(sigma_bg)
    eps_ambient = @(f) permittivity_water(f, sigma_bg(i));
    loglog(freqs, imag(eps_ambient(freqs)), '-.', 'LineWidth', lw, 'Color', color+[i/(lc/2) i/lc i/lc])
    hold on
end
grid on
% ylim([0 50])
xlabel("$f$ (Hz)",'Interpreter','latex','FontSize',fs)
ylabel("$\varepsilon_r''$",'Interpreter','latex','FontSize',fs)
% title("Permittivity of water (imaginary part)")
% hold off

eps_bfat = @(f) permittivity_tissue(2*pi*f, "Breast fat");
eps_btumor = @(f) permittivity_breastcancer(2*pi*f);
loglog(freqs, imag(eps_bfat(freqs)), 'LineWidth', lw, 'Color', [0.8500 0.3250 0.0980])
loglog(freqs, imag(eps_btumor(freqs)), 'LineWidth', lw, 'Color', [0.9290 0.6940 0.1250])
legend("$\sigma_w=$ "+sigma_bg(1),"$\sigma_w=$ "+sigma_bg(2),"$\sigma_w=$ "+sigma_bg(3),"$\sigma_w=$ "+sigma_bg(4), ...
    "Breast fat", "Breast cancer", 'interpreter','latex','fontsize',fs-2)
title("Permittivity (imaginary part)")
hold off


