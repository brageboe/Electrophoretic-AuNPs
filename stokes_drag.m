% References:
% [1] https://www.ias.ac.in/public/Volumes/pmsc/109/04/0441-0452.pdf

function beta = stokes_drag(eta, diameter, length)
    if diameter == length
        beta = stokes_drag_sphere(eta, diameter/2);
    elseif diameter > length
        beta = stokes_drag_oblate(eta, length, diameter);
    else
        beta = stokes_drag_prolate(eta, diameter, length);
    end
end

function beta = stokes_drag_sphere(eta, radius)
    beta_r = 6*pi*eta*radius;
    beta = [beta_r, beta_r, beta_r];
end

function beta = stokes_drag_oblate(eta, thickness, diameter)
% Equations (3.10) and (3.11) in [1]
    e = sqrt(1 - thickness^2/diameter^2);
    a = diameter/2;
    
    beta_ax = 8*pi*eta*a * e^3 / (e*sqrt(1-e^2) - (1 - 2*e^2)*asin(e));
    beta_tr = 16*pi*eta*a * e^3 / ((1 + 2*e^2)*asin(e) - e*sqrt(1-e^2));
    beta = [beta_tr, beta_tr, beta_ax];
end

function beta = stokes_drag_prolate(eta, diameter, length)
% Equations (3.5) and (3.7) in [1]
    e = sqrt(1 - diameter^2 / length^2);
    a = length/2;
    length = log((1+e)/(1-e));
    beta_ax = 16*pi*eta*a * e^3 / ((1 + e^2)*length - 2*e);
    beta_tr = 32*pi*eta*a * e^3 / (2*e + (3*e^2 - 1)*length);
    beta = [beta_tr, beta_tr, beta_ax];
end
