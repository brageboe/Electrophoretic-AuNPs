%% Permittivity due to electrophoretic movement of ellipsoidal particles suspended in a host medium
% omega: wave angular frequency [rad/s]
% phi: volume fraction of particles in suspension, phi = n * V_gnp /
% V_suspension
% eta: dynamic shear viscosity of host medium [N*s/m^2]
% epsilon_m: permittivity of host medium
% rho: density of particle material [kg/m^3]
% q: charge of a particle [C]
% D: diameter of particle [m]
% L: length of particle [m]
%
% Returns: epsilon = [epsilon_x, epsilon_y, epsilon_z], permittivity of
% suspension
%
% References:
% [1] https://www.ias.ac.in/public/Volumes/pmsc/109/04/0441-0452.pdf drag
% [2] DOI: 10.1088/0022-3727/45/7/075303, Sassaroli
% [3] DOI: 10.1088/1361-6463/aa7c8a, MD2017
%
% Written by Brage Bøe Svendsen, 2024
function epsilon_eph = permittivity_eph_suspension(freq, phi, eta, epsilon_host, q, D, L, options)
    arguments
       freq         % operating frequency
       phi          % volume fraction of gnps
       eta          % viscosity of host medium
       epsilon_host % permittivity of host matrix
       q            % net surface charge
       D            % diameter, gnp
       L            % length, gnp
       options.rho = 1.93e4; % mass density of metal core, default is Au
       options.includeAmbient = true; % false: purely electrophoretic (sassaroli)
       options.includeAmbientMD2017 = false; % true: MD2017 style
       options.printContrast = false;
    end
    
    beta = stokes_drag(eta, D, L);
    
    omega = 2*pi*freq;
    eps_0 = 8.8541878128e-12;
    V_GNP = (4/3) * pi * D * D * L / 8;
    m = options.rho * V_GNP;
     
    sigma_1 = phi * q^2 ./ (V_GNP * beta);
    tau_1 = m ./ beta;
    epsilon_flow = 1i/(eps_0*omega) * sigma_1./(1 - 1i*omega*tau_1);
    
    if options.includeAmbient && ~options.includeAmbientMD2017
        % Slightly modified from MD2017 eq (25): added the factor (1-phi).
        epsilon_eph = epsilon_flow + (1-phi)*epsilon_host;
    elseif options.includeAmbientMD2017
        % eq (25) from MD2017. 
        epsilon_eph = epsilon_flow + epsilon_host;
    else
        % Purely electrophoretic, Sassaroli eq (22).
        epsilon_eph = epsilon_flow;
    end
    if options.printContrast
        fprintf( "imag(eps_flow) - imag(eps_host) =\n\n\t" + (max(imag(epsilon_flow))-imag(epsilon_host)) + "\n")
    end
end

