function T = bioheat_uniform_water(freq, E0, t, suspension_radius, C_abs, host_permittivity, T0, res) 
    % Return temperature T in suspension of spherical volume.
    % Spatially uniform.
    % Ambient: water.
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
    c_w = 4186;     % Water specific heat capacity   [J/(kg*K)] 
    rho_w = 998;    % Water density                  [kg/m^3]   
    k_w = 0.6071;	% Water thermal conductivity     [W/(m*K)]  
    
    omega = 2*pi*freq;
    k0 = omega/299792458;
    mu0 = 1.25663706127e-6;     % vacuum permeability
    suspension_volume = (4/3)*pi*suspension_radius^3;

    %res = 3; % minimum value for grid; we only want the max value.
    T0_grid = T0.*ones(res, res, res);	% Initial temperature (3d matrix) [degC]
    S = ones(res, res, res);            % Heatrate due to external heat source (RF)

    k_host = k0*sqrt(host_permittivity);     
    I0 = (E0.*E0)*real(k_host)/(2*omega*mu0);
    
    w_abs = C_abs .* I0 / suspension_volume; 
    S = S .* w_abs/(rho_w*c_w); 
    
    region_side_length = suspension_volume ^ (1/3) * 0.01; 
    dx = region_side_length / res; 
    
    D = k_w/(rho_w*c_w);
    P = 0;  % blood perfusion rate
    Ta = 0; % arterial blood temperature

    T = bioheatExact(T0_grid, S, [D, P, Ta], dx, t);  
end