% Sassaroli: https://iopscience.iop.org/article/10.1088/0022-3727/45/7/075303/pdf
clear variables 

e0 = 1.6e-19; % Elementary charge
eta = 1e-3; % Viscosity of host medium
rho = 1.93e4; % Particle material density

% R = 5e-9;
phi = 2e-3; % Volume fraction of particles:tissue
% D = 10e-9;
% L = 10e-9;

f = 1e9;

Q = @(R) e0 * (3 * 1e9*R + 0.5 * (1e9*R)^2); % Sassaroli expression for excess charge
q = @(D, L) 0.2*surface_area(D, L); % Constant surface charge density
q2 = @(D, L) Q((D*D*L)^0.3333 / 2); % Sassaroli expression using "effective radius"
q3 = @(D, L) 500*e0; % Constant charge

epsilon_m = @(omega) permittivity_tissue(omega, "Breast fat");
epsilon_eph = @(f, phi, D, L) permittivity_eph_suspension(f, phi, eta, epsilon_m(f*2*pi), q2(D,L), D, L);
% epsilon_eph_1 = @(f, phi, R) permittivity_eph_suspension(f, phi, eta, 1, rho, Q(R), R*2, R*2);

R_s = 1e-3; % Arbitrary
C_abs = @(f, phi, D, L) mean(C_abs_spheroid(f, epsilon_eph(f, phi, D, L), epsilon_m(f*2*pi), R_s*2, R_s*2, approximation_check=false));
C_amb = @(f) C_amb_spheroid(f, epsilon_m(f*2*pi), R_s*2, R_s*2);
F_abs = @(f, phi, D, L) C_abs(f, phi, D, L) / C_amb(f);

% Sassaroli plot
%figure();
%set_standard_plot_style();
%fplot(@(R) imag(mean(epsilon_eph_1(15e6, 2e-2, R))), [2e-9, 15e-9]); 
%xlabel("R (m)");
%ylabel("\epsilon''");

figure();
make_plot_2d(@(D, L) F_abs(f, phi, D, L), [1e-9 1e-4], [1e-9 1e-4], zlabel="$F_{abs}$", res=200, y_log=true);
set_standard_plot_style();
xlabel("$D$ (m)");
ylabel("$L$ (m)");

figure();
make_plot_2d(@(D,L) q(D,L)/e0, [1e-9 1e-4], [1e-9 1e-4], zlabel="Net surface charge ($e_0$)");
set_standard_plot_style();
xlabel("$D$ (m)");
ylabel("$L$ (m)");

function A = surface_area(D, L)
    if D == L
        A = pi*D^2;
    else
        e = sqrt(1 - D^2/L^2);
        A = 0.5*pi*D^2 + (0.5*pi*D*L/e)*asin(e);
    end
end