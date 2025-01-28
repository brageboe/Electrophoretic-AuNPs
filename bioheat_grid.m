%% TO-DOs
% - does not agree with other codes e.g. example_spheroidalGNPs_DL.m
%   - almost zero temperature increase even after 5 min
% - temperature (incl. center of suspension) changes wrt volumes

function T = bioheat_grid(frequency, E0, exposure_time, suspension_volume, ambient_volume, suspension_permittivity, ambient_permittivity, temperature_0, res, tissue)
% Returns a [res x res x res] cubical grid of temperature T after exposure_time seconds
    arguments
        frequency
        E0  % electric field strength [V/m]
        exposure_time
        suspension_volume
        ambient_volume
        suspension_permittivity
        ambient_permittivity
        temperature_0
        res
        tissue = "Breast tumor" % ambient and host media must be same 
    end
    % Tissue data from Table 8 in [1] Camilleri et al. DOI: 10.3390/s22103894
    if tissue == "Breast fat"
        rho_t = 932;        % mass density                      [kg/m3]
        c_t = 2220;         % specific heat capacity            [J/kg/K]
        k_eff = 0.171;      % effective thermal conductivity    [W/m/K]
        omega_t = 0.886;    % blood perfusion rate              [kg/s/m3] 
    elseif tissue == "Breast tumor"
        rho_t = 1066;       % mass density                      [kg/m3]
        c_t = 3610;         % specific heat capacity            [J/kg/K]
        k_eff = 0.511;      % effective thermal conductivity    [W/m/K]
        omega_t = 4.830;    % blood perfusion rate              [kg/s/m3]
    elseif tissue == "Water"
        rho_t = 998;        % Water density                     [kg/m^3] 
        c_t = 4186;         % Water specific heat capacity      [J/(kg*K)] 
        k_eff = 0.6071;     % 0.6071 @ 25degC; 0.62856 @ 40degC.[W/m/K]
        omega_t = 0;
    else
        error("Input 'tissue' may only take value 'Breast fat', 'Breast tumor', or 'Water'.")
    end
    c_b = 3622.50;          % specific heat capacity of blood [J/kg/K]. Table 2 [1].
    % rho_b = 1049.75;      % mass density [kg/m3]. not needed; omega_t=rho_b*w_t

    omega = 2*pi*frequency;     % vacuum angular frequency
    k0 = omega/299792458;       % vacuum wavenumber
    mu0 = 1.25663706127e-6;     % vacuum permeability
    % suspension_volume = (4/3)*pi*suspension_radius^3;
    suspension_radius = (3*suspension_volume/(4*pi))^(1/3);

    T0_grid = temperature_0.*ones(res, res, res);	% Initial temperature (3d matrix) [degC]
    rf_power_absorption = ones(res, res, res);      % Heatrate due to external heat source [W/m3]

    k_amb = k0*sqrt(ambient_permittivity);          % ambient wavenumber
    I0 = (E0*E0) * real(k_amb) / (2*omega*mu0);     % intensity of incident wave  
    
    C_abs = mean(C_abs_spheroid(frequency, suspension_permittivity, ambient_permittivity, suspension_radius*2, suspension_radius*2, approximation_check=true, report_warning=true));
    C_amb = C_amb_spheroid(frequency, ambient_permittivity, suspension_radius*2, suspension_radius*2);

    w_amb = C_amb * I0 / suspension_volume;
    w_abs = C_abs * I0 / suspension_volume; 
    
    domain_side_length = ambient_volume ^ (1/3);    % side length of cube shaped ambient domain [m]
    dx = domain_side_length / res;                  % distance between grid points
    
    % Loop through every point on the grid and fill in the power absorption 
    % per unit volume. If the point belongs to the bolus injection use
    % w_np, otherwise use w_amb
    for i = 1:res
        for j = 1:res
            for k = 1:res
                pos = ([i, j, k] - 0.5) * dx - domain_side_length / 2; % convert matrix index to actual point in space, origo in center of the cube
                if (norm(pos) <= suspension_radius)
                    rf_power_absorption(i, j, k) = w_abs;
                else
                    rf_power_absorption(i, j, k) = w_amb;
                end
            end
        end
    end
    
    % as in http://www.k-wave.org/documentation/bioheatExact.php
    D = k_eff / (rho_t * c_t);          % spatial heat conduction factor
    P = c_b * omega_t / (rho_t * c_t);  % blood perfusion rate
    Ta = 37;                            % arterial blood temperature
    S = rf_power_absorption / (rho_t * c_t);

    T = bioheatExact(T0_grid, S, [D, P, Ta], dx, exposure_time);  
end