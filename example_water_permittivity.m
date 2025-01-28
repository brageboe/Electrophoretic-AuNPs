% Permittivity (imaginary part) of water for different ionic conductions.

sigma_bg = [0,1e-3,1e-2,1e-1];       
freqs = logspace(5, 10, 1e3); 

figure();
lw = 1.5;
fs = 12;
for i=1:length(sigma_bg)
    eps_ambient = @(f) permittivity_water(f, sigma_bg(i));
    semilogx(freqs, imag(eps_ambient(freqs)), 'LineWidth', lw)
    hold on
end
grid on
ylim([0 50])
xlabel("$f$ (Hz)",'Interpreter','latex','FontSize',fs)
ylabel("$\varepsilon_r''$",'Interpreter','latex','FontSize',fs)
legend("$\sigma_w=$ "+sigma_bg(1),"$\sigma_w=$ "+sigma_bg(2),"$\sigma_w=$ "+sigma_bg(3),"$\sigma_w=$ "+sigma_bg(4),'interpreter','latex','fontsize',fs-2)
title("Permittivity of water (imaginary part)")
hold off