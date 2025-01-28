% sassoroli_fig8
clear variables
clear figures

e0 = 1.6e-19; % Elementary charge
eta = 1e-3; % Viscosity of host medium
% rho = 1.93e4; % Particle material density


q = @(R) e0*(3*1e9*R + 0.5*(1e9*R).^2); % Net surface charge, Sassaroli eq (25)
% epsilon_m = @(f) permittivity_water(f);
epsilon_eph = @(f, R, phi, Q) permittivity_eph_suspension(f, phi, eta, 1, Q, R*2, R*2, includeAmbient=false);

figure()
%% Fig 8a: wrt frequency, f
freq = linspace(1,100)*1e6;
a = 5e-9;
N = [27.5 20]; % note that here q(a)=27.5e0
phi = 2e-2;     

eps_flow_f1 = zeros(3, length(freq));
eps_flow_f2 = zeros(3, length(freq));
for i=1:length(freq)
    eps_flow_f1(:,i) = epsilon_eph(freq(i), a, phi, N(1)*e0);%permittivity_eph_suspension(freq(i), phi, eta, epsilon_m(freq), rho, N(1)*e0, a*2, a*2, includeAmbient=false);%epsilon_eph(freq(i), a, phi);
    eps_flow_f2(:,i) = epsilon_eph(freq(i), a, phi, N(2)*e0);%permittivity_eph_suspension(freq(i), phi, eta, epsilon_m(freq), rho, N(2)*e0, a*2, a*2, includeAmbient=false);
end

% figure()
subplot(1,3,1)
plot(freq/1e6, imag(eps_flow_f1(1,:)))
hold on
plot(freq/1e6, imag(eps_flow_f2(1,:)))
legend("$N=$"+N(1), "$N=$"+N(2), 'interpreter', 'latex')
xlabel("Frequency [MHz]", 'interpreter', 'latex')
ylabel("$\varepsilon_{flow}''$", 'interpreter', 'latex')
hline = findobj(gcf, 'type', 'line');
set(hline, 'LineWidth', 1.5)
grid on 
hold off

%% Fig 8b: wrt volume fraction, phi
freq = 15e6;
a = [5e-9 25e-9];
phi = logspace(-6, -1);     

eps_flow_p1 = zeros(3, length(phi));
eps_flow_p2 = zeros(3, length(phi));
for i=1:length(phi)
    eps_flow_p1(:,i) = epsilon_eph(freq, a(1), phi(i), q(a(1)));
    eps_flow_p2(:,i) = epsilon_eph(freq, a(2), phi(i), q(a(2)));
end

% figure()
subplot(1,3,2)
semilogx(phi, imag(eps_flow_p1(1,:)))
hold on
semilogx(phi, imag(eps_flow_p2(1,:)))
legend("$a=$"+a(1)/1e-9+" nm", "$a=$"+a(2)/1e-9+" nm", 'interpreter', 'latex')
xlabel("$\phi$", 'interpreter', 'latex')
ylabel("$\varepsilon_{flow}''$", 'interpreter', 'latex')
xlim([min(phi) 2e-2])
hline = findobj(gcf, 'type', 'line');
set(hline, 'LineWidth', 1.5)
grid on 
hold off

%% Fig 8c: wrt gnp radius, a
freq = 15e6; 
a = linspace(1,15)*1e-9;
phi = [2e-2, 2e-3]; 

eps_flow_a1 = zeros(3, length(a));
eps_flow_a2 = zeros(3, length(a));
charge = q(a);
for i=1:length(a)
    eps_flow_a1(:,i) = epsilon_eph(freq, a(i), phi(1), charge(i));
    eps_flow_a2(:,i) = epsilon_eph(freq, a(i), phi(2), charge(i));
end
clear charge

% figure(3)
subplot(1,3,3)
hold on
plot(a/1e-9, imag(eps_flow_a1(1,:)))
plot(a/1e-9, imag(eps_flow_a2(1,:)))
legend("$\phi=$"+phi(1), "$\phi=$"+phi(2), 'interpreter', 'latex')
xlabel("Radius [nm]", 'interpreter', 'latex')
ylabel("$\varepsilon_{flow}''$", 'interpreter', 'latex')
hline = findobj(gcf, 'type', 'line');
set(hline, 'LineWidth', 1.5)
grid on 
hold off
