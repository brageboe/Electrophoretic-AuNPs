function T = bioheat_uniform_breasttumor(freq, E0, t, suspension_radius, C_abs, host_permittivity, T0, res) 
    arguments
        freq
        E0
        t
        suspension_radius
        C_abs
        host_permittivity
        T0
        res = 3
    end
    % Tissue data from Table 8 in [1] Camilleri et al. DOI: 10.3390/s22103894
    rho_t = 1066;       % mass density                      [kg/m3]
    c_t = 3610;         % specific heat capacity            [J/kg/K]
    k_eff = 0.511;      % effective thermal conductivity    [W/m/K]
    omega_t = 4.830;    % blood perfusion rate              [kg/s/m3]
    c_b = 3622.50;      % specific heat capacity of blood [J/kg/K]. Table 2 [1].

    omega = 2*pi*freq;
    k0 = omega/299792458;
    mu0 = 1.25663706127e-6;     % vacuum permeability
    suspension_volume = (4/3)*pi*suspension_radius^3;

    % res = 3; % minimum value for grid; we only want the max value.
    T0_grid = T0.*ones(res, res, res);	% Initial temperature (3d matrix) [degC]
    S = ones(res, res, res);            % Heatrate due to external heat source (RF)
    
    k_host = k0*sqrt(host_permittivity);     
    I0 = (E0.*E0)*real(k_host)/(2*omega*mu0);

    w_abs = C_abs .* I0 / suspension_volume; 
    S = S .* w_abs/(rho_t*c_t); 

    region_side_length = suspension_volume ^ (1/3) * 0.01; 
    dx = region_side_length / res; 

    % as in http://www.k-wave.org/documentation/bioheatExact.php
    D = k_eff / (rho_t * c_t);          % spatial heat conduction factor
    P = c_b * omega_t / (rho_t * c_t);  % blood perfusion rate
    Ta = 37;                            % arterial blood temperature

    T = bioheatExact(T0_grid, S, [D, P, Ta], dx, t);  
end