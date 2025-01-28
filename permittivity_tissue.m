% Returns function of omega, basically koden som vi fick av brage
function epsilon = permittivity_tissue(omega, tissue)
    %% 4-term-Cole-Cole function for various human tissues
    % Input table should contain data from http://niremf.ifac.cnr.it/docs/DIELECTRIC/AppendixC.html#FF
    
    %%
    % Inputs: tissue (e.g. 'Breast fat') and frequency spectrum omega.
    persistent tissue_name table epsilon_func;
    %class(tissue_name)
    if isempty(tissue_name) || tissue_name ~= tissue
        tissue_name = tissue;
        if isempty(table)
            table = importdata("Tissue Data/summaryTable_colecole.txt"); % Imports the data in summaryTable
        end
        index = find(strcmp(table.textdata, tissue_name)); % Finds the data in summaryTable related to the choosen tissue
        tissueData = table.data(index, :); % Stores the parameters for the Cole-Cole equation 

        epsinf = tissueData(1); % High water content
        Deltaeps = [tissueData(2); tissueData(5); tissueData(9); tissueData(12)]; % Magnitude of dispersion
        tau = [tissueData(3)*(10^(-12)); tissueData(6)*(10^(-9)); tissueData(10)*(10^(-6)); tissueData(13)*(10^(-3))]; % Relaxation time
        alpha = [tissueData(4); tissueData(7); tissueData(11); tissueData(14)]; % Distribution parameter
        sigma = tissueData(8); % Static ionic conductivity
        eps_0 = 8.8541878128e-12;% vacuum permittivity
        cond = @(omega) sigma./(1i*eps_0*omega); % Conductivity term
        
        
        %eps_disp = Deltaeps./(1+1i*omega.*tau).^(1-alpha); % Permittivity of each separate Debye dispersion region
        eps = @(omega) sum(Deltaeps./(1+(1i*omega.*tau).^(1-alpha)));
        
        % Complex relative permittivity for surrounding medium with cole-cole equation:
        epsilon_func = @(omega) conj(epsinf + eps(omega) + cond(omega));
    end

    epsilon = epsilon_func(omega);
end