% D, L dimensions of spheroid
% Symmetric around z-axis; L1 = L2
function L_i = geometric_factors_spheroid(D, L)
    if D == L
        L_i = [1, 1, 1] / 3; %Bohren & Huffman, p 146
    elseif (D > L)
        L_i = geometric_factors_oblate_spheroid(D, L);
    else
        L_i = geometric_factors_prolate_spheroid(D, L);
    end
end
%(5.33) in Bohren & Huffman
function L_i = geometric_factors_prolate_spheroid(D, L)
    e = sqrt(1 - D^2 / L^2);
    L3 = ((1 - e^2)/e^2) * (-1 + (1/(2*e))*log((1 + e) / (1 - e))); 
    L1 = (1 - L3) / 2; %Bohren & Huffman, p 146: L1 + L2 + L3 = 1
    L_i = [L1, L1, L3];
end
% (5.34) in Bohren & Huffman
function L_i = geometric_factors_oblate_spheroid(D, T)
    e = sqrt(1 - T^2 / D^2);
    g = sqrt((1 - e^2) / e^2);
    L1 = (g / (2*e^2))*(pi/2 - atan(g)) - g^2 / 2;
    L3 = 1 - 2*L1; %Bohren & Huffman, p 146: L1 + L2 + L3 = 1
    L_i = [L1, L1, L3];
end