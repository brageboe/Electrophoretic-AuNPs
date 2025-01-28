function q = net_charge_Rostalski(D,L)
% Rostalski J and Quinten M 1996 Colloid Polym. Sci. 274 648
% (found in Sassaroli2012 as ref [51])
    arguments
        D   % diameter, gnp
        L   % length, gnp
    end
    e0 = 1.602176634e-19;      % electron charge
    Q = @(R) e0 * (3 * 1e9*R + 0.5 * (1e9*R)^2); % Sassaroli eq (25) for excess charge
    if D == L % the two ifs are the same but first is more accurate for spheres
        q = Q(D/2);
    else
        q = Q((D*D*L)^0.33333 / 2); % "effective spherical radius" for spheroids
    end
end