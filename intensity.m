function I = intensity(freq, E0, permittivity)
    omega = 2*pi*freq;
    k0 = omega/299792458;
    mu0 = 1.25663706127e-6;
    k_host = k0 * sqrt(permittivity);
    I = (E0.*E0)*real(k_host)/(2*omega*mu0);
end
