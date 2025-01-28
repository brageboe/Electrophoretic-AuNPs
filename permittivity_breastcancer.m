function epsilon = permittivity_breastcancer(omega)
    % CANCER PERMITTIVITY 
    % with Cole parameters (from Lazebnik)
    % Available: https://iopscience.iop.org/article/10.1088/0031-9155/52/20/002/pdf
    
    eps_0 = 8.8541878128e-12;   % vacuum permittivity
    
    sigma_C = 0.794;
    Deltaeps_C = 50.09;
    tau_C = 10.50e-12;
    alpha_C = 0.051;
    epsinf_C = 6.749;

    cond_C = sigma_C./(1i*eps_0*omega); % Conductivity term
    eps_disp_C = Deltaeps_C./(1+(1i*omega.*tau_C).^(1-alpha_C)); % Permittivity of each separate Debye dispersion region
    epsilon = conj(epsinf_C + eps_disp_C + cond_C); % Complex relative permittivity
end