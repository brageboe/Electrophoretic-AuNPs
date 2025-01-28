function q = net_charge_Gupta(D,L,n)
% GNP surface charge as determined by a constant [charge per surface area].
% Inspiried by Gupta and Rai 2016 (DOI: 10.1038/srep45292), where n=5.
    arguments 
        D       % diameter, gnp (m)
        L       % length, gnp (m)
        n = 5   % number of elementary charges per nm surface area (#/nm^2)
    end
    e0 = 1.602176634e-19;      % electron charge
    q = e0 * surface_area(D,L) * n * 1e18;
end


function A = surface_area(D, L)
    if D == L
        A = pi*D.^2;
    else
        e = sqrt(1 - D^2/L^2);
        A = 0.5*pi*D^2 + (0.5*pi*D*L/e)*asin(e);
    % elseif L < D % oblate
    %     e = sqrt(1 - L^2/D^2);
    %     A = 0.5*pi*D^2 + (0.5*pi*D^2)*atanh(e)*(1-e^2)/e;
    % Expression for oblate can be used for prolates, and vice versa.
    % Difference is e becomes imaginary and is no longer the eccentricity.
    end
end