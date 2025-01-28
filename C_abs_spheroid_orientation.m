function C_abs = C_abs_spheroid_orientation(freq, epsilon_p, epsilon_m, D_p, L_p, options)
% Similar to C_abs_spheroid.m but here you have options on GNP orientations.
% flow_oriented: returns C_abs for which friction resistance is lowest.
% field_oriented: returns C_abs for which suspension polarizability is highest.
% default: returns C_abs = [C_abs_1 C_abs_2 C_abs_3] as normal.
%       -> randomly oriented by then taking the average mean(C_abs)
% Written by Brage Bøe Svendsen, September 2024.
%
% Todo:
% For cases with max(alpha) or min(beta) in multiple directions, instead 
% of arbitrarily choosing one, check whether one of them has smaller beta
% or larger alpha, respectively, and choose that. Since it is conceivable
% that that direction will be favourable for the GNPs to orient along.
% However, in the cases I have compared, the direction with max(alpha) is
% also the direction of min(beta), so this might be relevant only for
% obscure cases.
    arguments
        freq            % frequency
        epsilon_p       % scattering particle
        epsilon_m       % surrounding medium
        D_p             % scattering particle diameter
        L_p             % scattering particle length
        options.approximation_check = true
        options.report_warning = false
        options.field_oriented = false
        options.flow_oriented = false 
        options.D_GNP
        options.L_GNP
        options.eta
    end
    omega = 2 * pi * freq;
    c0 = 299792458; 
    k0 = omega / c0;
    kp = k0 * sqrt(epsilon_p);
    km = k0 * sqrt(epsilon_m);
    
    a = [D_p, D_p, L_p] / 2;
    L_i = geometric_factors_spheroid(D_p, L_p);
    V = (4/3)*pi*a(1)*a(2)*a(3);
    
    if options.field_oriented
        alpha = polarizability(epsilon_p, epsilon_m, V, L_i);
        direction_i = find(alpha == max(alpha));
        % cases where alpha has maximum in multiple directions 
        % (oblates and spheres)
        if length(direction_i) > 1    
            direction_i = direction_i(1);
        end
        C_abs = C_abs_direction(k0, epsilon_p, epsilon_m, V, alpha, direction_i);
    elseif options.flow_oriented
        alpha = polarizability(epsilon_p, epsilon_m, V, L_i);
        beta = stokes_drag(options.eta, options.D_GNP, options.L_GNP);
        direction_i = find(beta == min(beta));
        % cases where beta has minimum in multiple directions 
        % (oblates and spheres)
        if length(direction_i) > 1
            direction_i = direction_i(1);
        end
        C_abs = C_abs_direction(k0, epsilon_p, epsilon_m, V, alpha, direction_i);
    else
        alpha = polarizability(epsilon_p, epsilon_m, V, L_i);
        numerator = k0 * abs(epsilon_m).^2 .* imag(epsilon_p) .* abs(alpha).^2;
        denominator = V * real(sqrt(epsilon_m)) .* abs(epsilon_p - epsilon_m).^2;
        C_abs_i =  numerator ./ denominator;
    end

    field_to_wave_direction = [2, 3, 1]; % See Rommelfanger2021, Fig1b. S in direction 1 => E in direction 2, etc.
    small_threshold = 0.1;

    if options.field_oriented || options.flow_oriented
        check = [real(kp), imag(kp), real(km), imag(km)]*a(field_to_wave_direction(direction_i)); % Rommelfanger2021, section IV
        if max(check) > small_threshold && options.approximation_check
            C_abs = 0;
            if options.report_warning == true
                warning("Bounds for electrostatic approximation exceeded @ D="+D_p+", L="+L_p+", f="+omega/(2*pi)+", axis i="+direction_i+".")
            end
        end
    else
        C_abs = [0, 0, 0];
        for i=1:3
            check = [real(kp), imag(kp), real(km), imag(km)]*a(field_to_wave_direction(i)); % Rommelfanger2021, section IV
            if max(check) <= small_threshold || ~options.approximation_check
               C_abs(i) = C_abs_i(i);
            elseif options.report_warning == true 
                if max(check) > small_threshold && options.approximation_check == true 
                    warning("Bounds for electrostatic approximation exceeded @ D="+D_p+", L="+L_p+", f="+omega/(2*pi)+", axis i="+i+".")
                end
            end
        end
    end
end
%%
function C_abs = C_abs_direction(k0, epsilon_p, epsilon_m, volume, polarizability, direction)
% Sets C_abs equal to C_abs in the given direction.
    numerator = k0 * abs(epsilon_m).^2 .* imag(epsilon_p(direction)) .* abs(polarizability(direction)).^2;
    denominator = volume * real(sqrt(epsilon_m)) .* abs(epsilon_p(direction) - epsilon_m).^2;
    C_abs = max(numerator ./ denominator); 
end

function alpha = polarizability(epsilon_p, epsilon_m, volume, geometric_factors)
% geometric_factors = [L_1 L_2 L_3].
% epsilon_p = [eps_1 eps_2 eps_3] or eps_scalar.
    numerator = 3 * volume * (epsilon_p - epsilon_m);
    denominator = 3*epsilon_m + 3*(epsilon_p - epsilon_m).*geometric_factors;
    alpha = numerator ./ denominator;
end