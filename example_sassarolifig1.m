%% example_sassarolifig1
f = linspace(1,100,100).*1e6;
sigma = [0 8e-4 16e-4 32e-4];

eps_w = @(f, sigma) permittivity_water(2*pi*f, sigma);

figure(1)
plot(f/1e6,imag(eps_w(f, sigma(1))))
hold on
plot(f/1e6,imag(eps_w(f, sigma(2))))
plot(f/1e6,imag(eps_w(f, sigma(3))))
plot(f/1e6,imag(eps_w(f, sigma(4))))

legend("$\sigma=$"+sigma(1)+" S/m","$\sigma=$"+sigma(2)+" S/m", ...
    "$\sigma=$"+sigma(3)+" S/m","$\sigma=$"+sigma(4)+" S/m", 'interpreter', 'latex')
xlabel("Frequency [MHz]", 'interpreter', 'latex')
ylabel("$\varepsilon_{water}''$", 'interpreter', 'latex')
hline = findobj(gcf, 'type', 'line');
set(hline, 'LineWidth', 1.5)
grid on 
hold off

%% over large spectrum
f = logspace(3,10,100); % 1 kHz to 10 GHz

figure(2)
loglog(f,imag(eps_w(f, sigma(1))))
hold on
loglog(f,imag(eps_w(f, sigma(2))))
loglog(f,imag(eps_w(f, sigma(3))))
loglog(f,imag(eps_w(f, sigma(4))))

legend("$\sigma=$"+sigma(1)+" S/m","$\sigma=$"+sigma(2)+" S/m", ...
     "$\sigma=$"+sigma(3)+" S/m","$\sigma=$"+sigma(4)+" S/m", 'interpreter', 'latex')
xlabel("Frequency [MHz]", 'interpreter', 'latex')
ylabel("$\varepsilon_{water}''$", 'interpreter', 'latex')
hline = findobj(gcf, 'type', 'line');
set(hline, 'LineWidth', 1.5)
grid on 