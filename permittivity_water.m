function epsilon = permittivity_water(freq, sigma)
% Debye permittivity of water as given in Sassaroli, but instead with exp(-iwt) convention.
    arguments
        freq
        sigma = 0; % Ionic contribution to total conductivity. Zero for pure water.
    end
    omega = 2*pi*freq;
    eps_0 = 8.8541878128e-12;
    eps_inf = 4.5; % water upper limit
    eps_static = 78.3; % water lower limit
    Delta_eps = eps_static - eps_inf;
    f_c = 19.5e9; % water characteristic frequency
    tau = 1 / (2*pi*f_c); % relaxation time
    
    epsilon = eps_inf + Delta_eps./(1-1i*omega*tau) - sigma./(1i*omega*eps_0);
end