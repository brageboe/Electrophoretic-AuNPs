function C_abs = C_abs_spheroid(freq, epsilon_p, epsilon_m, D, L, options)
% Written by Olle Hennert 2023. Comments and edits by Brage Bøe Svendsen.
    arguments
        freq        % frequency
        epsilon_p   % scattering particle
        epsilon_m   % surrounding medium
        D           % scattering particle diameter
        L           % scattering particle length
        options.approximation_check = true
        options.report_warning = false
    end
    omega = 2 * pi * freq;
    c0 = 299792458; 
    k0 = omega / c0;
    kp = k0 * sqrt(epsilon_p);
    km = k0 * sqrt(epsilon_m);
    
    a = [D, D, L] / 2;
    L_i = geometric_factors_spheroid(D, L);
    V = (4/3)*pi*a(1)*a(2)*a(3);
    numerator = k0*V*abs(epsilon_m).^2.*imag(epsilon_p);
    denominator = L_i.^2.*real(sqrt(epsilon_m)).*abs(epsilon_p - epsilon_m.*(L_i-1)./L_i).^2;
    C_abs_i =  numerator ./ denominator;

    field_to_wave_direction = [2, 3, 1]; % See Rommelfanger2021, Fig1b. S in direction 1 => E in direction 2, etc.
    small_threshold = 0.1;
    C_abs = [0, 0, 0];
    for i=1:3
        check = [real(kp), imag(kp), real(km), imag(km)]*a(field_to_wave_direction(i)); % Rommelfanger2021, section IV
        if max(check) <= small_threshold || options.approximation_check == false
           C_abs(i) = C_abs_i(i);
        elseif options.report_warning == true % added 2/8/24, BBS
            if max(check) > small_threshold && options.approximation_check == true % added 4/6/2024, BBS
                warning("Bounds for electrostatic approximation exceeded @ D="+D+", L="+L+", f="+omega/(2*pi)+", axis i="+i+".")
            end
        end
    end
end