clear variables

%% Parameters
e0 = 1.602176634e-19;     
omega = 2*pi*13.56e6;
a = 1.65e-9;                    % GNP radius
c_w = 4186;                     % Water specific heat capacity   [J/(kg*K)] 
rho_w = 998;                    % Water density                  [kg/m^3]  
rho_Au = 19300;                 % Gold density
eta = 0.87e-3;                  % Water viscosity
V_s = 2e-6;                     % Sample volume, 2 mL
V_p = (4/3)*pi*a^3;             % GNP volume
beta = 6*pi*eta*a;              % Friction constant
tau = 2*rho_Au*a*a/(9*eta);     % Relaxation time

heatrate = [0.0218, 0.0277, 0.0311];    % using experimental (background-subtracted) heatrates averaged from example_Collins2018_FigS8.m
q = [-5.9*e0, -6.7*e0, -7.4*e0];        % net charge per GNP for samples [A,B,C]

%% Function handle & function space
phifun = @(E, sample) 2*heatrate(sample)*c_w*rho_w*V_p*beta*(1+tau^2*omega^2)./(q(sample)^2*E.^2);
% Efun = @(phi, sample) sqrt(2*heatrate(sample)*c_w*rho_w*V_p*beta*(1+tau^2*omega^2)./(q(sample)^2*phi));

E = linspace(1e0,1e4,1000);
% phi = logspace(-7,-1,1000);

%% Figures
figure()
lw = 1.5;
semilogx(phifun(E,1), E, '-', 'LineWidth',lw)
hold on
semilogx(phifun(E,2), E, '--', 'LineWidth',lw)
semilogx(phifun(E,3), E, '-.', 'LineWidth',lw)
legend("A","B","C")
xlabel("$\phi$", 'interpreter', 'latex')
ylabel("$E$ [V/m]", 'interpreter', 'latex')
xlim([8e-4 1e-1])
% ylim([1 1e3])
grid on
hold off


% figure()
% loglog(phi, Efun(phi,1), '-', 'LineWidth',lw)
% hold on
% loglog(phi, Efun(phi,2), '--', 'LineWidth',lw)
% loglog(phi, Efun(phi,3), '-.', 'LineWidth',lw)
% legend("A","B","C")
% ylabel("$E$ [V/m]", 'interpreter', 'latex')
% xlabel("$\phi$", 'interpreter', 'latex')
% xlim([1e-7 1e-1])
% ylim([1 1e3])
% grid on
% hold off