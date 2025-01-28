function C_amb = C_amb_spheroid(freq, epsilon_m, D, L)
    c0 = 299792458; %speed of light in vacuum
    omega = 2 * pi * freq;
    k0 = omega / c0;
    a = [D, D, L] / 2;
    V = (4/3)*pi*a(1)*a(2)*a(3);
    C_amb = 2*k0*V*imag(sqrt(epsilon_m));
end