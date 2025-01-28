clear variables 
clc
A = linspace(3, 200, 1e3).*1e-9;  % length of one gnp axis (diameter or length)
B = min(A); % prolates B=diameter; oblates B=thickness; spheres B=A=diameter
eta = 0.87e-3;

% SPHERES
beta_sphere = zeros(3,length(A));
for i=1:length(A)
    beta_sphere(:,i) = stokes_drag(eta, A(i), A(i));
end
% PROLATES
beta_prolate = zeros(3,length(A));
for i=1:length(A)
    beta_prolate(:,i) = stokes_drag(eta, B, A(i));
end
% OBLATES
beta_oblate = zeros(3,length(A));
for i=1:length(A)
    beta_oblate(:,i) = stokes_drag(eta, A(i), B);
end

ymin = min([mean(beta_sphere) mean(beta_prolate) mean(beta_oblate)]);
ymax = max([mean(beta_sphere) mean(beta_prolate) mean(beta_oblate)]);
%% FIGURE: average beta
set_standard_plot_style();
f1 = figure();
subplot(1,3,1)
semilogy(A/1e-9, mean(beta_sphere))
title("Spheres")
xlabel("Diameter [nm]")
% ylabel("$\frac{1}{3}\sum_i^3 \beta_i$ (mPa$\cdot$s)")
ylabel("mean($\beta_i$) [Pa$\cdot$s$\cdot$m]")
ylim([ymin ymax])
grid on
% set_standard_plot_style();
% hold off

subplot(1,3,2)
semilogy(A/1e-9, mean(beta_prolate))
title("Prolates")
xlabel("Length [nm]")
ylim([ymin ymax])
grid on
% set_standard_plot_style();
% hold off

subplot(1,3,3)
semilogy(A/1e-9, mean(beta_oblate))
title("Oblates")
xlabel("Diameter [nm]")
ylim([ymin ymax])
grid on
% set_standard_plot_style();
% hold off

axesHandles = get(gcf, 'children');
set(axesHandles, 'LineWidth', 1);
set(axesHandles, 'FontSize', 12);
set(axesHandles(1), 'YTickLabel', []);
set(axesHandles(2), 'YTickLabel', []);

f1.Position = [680 458 560 280]; % [x0 y0 width height]

%% FIGURE: axial beta
ymin = min([min(min(beta_sphere)) min(min(beta_sphere)) min(min(beta_sphere))]);
ymax = max([max(max(beta_sphere)) max(max(beta_sphere)) max(max(beta_sphere))]);

f2 = figure();
subplot(1,3,1)
semilogy(A/1e-9, beta_sphere(1,:))
title("Spheres")
xlabel("Diameter [nm]")
% ylabel("$\frac{1}{3}\sum_i^3 \beta_i$ (mPa$\cdot$s)")
ylabel("$\beta_i$ [Pa$\cdot$s$\cdot$m]")
ylim([ymin ymax])
grid on
% legend("$\beta_r$", 'location', 'southeast')

subplot(1,3,2)
semilogy(A/1e-9, beta_prolate(1,:))
hold on
% semilogy(A/1e-9, beta_prolate(2,:))
semilogy(A/1e-9, beta_prolate(3,:))
title("Prolates")
xlabel("Length [nm]")
ylim([ymin ymax])
grid on
% legend("$\beta_{tr}$", "$\beta_{ax}$", 'location', 'southeast')
hold off

subplot(1,3,3)
semilogy(A/1e-9, beta_oblate(1,:))
hold on
% semilogy(A/1e-9, beta_oblate(2,:))
semilogy(A/1e-9, beta_oblate(3,:))
title("Oblates")
xlabel("Diameter [nm]")
ylim([ymin ymax])
grid on
% legend("$\beta_{tr}$", "$\beta_{ax}$", 'location', 'southeast')
hold off

axesHandles = get(gcf, 'children');
set(axesHandles, 'LineWidth', 1);
set(axesHandles, 'FontSize', 12);
set(axesHandles(1), 'YTickLabel', []);
set(axesHandles(2), 'YTickLabel', []);
legend("$\beta_{tr}$", "$\beta_{ax}$", 'location', 'southeast')

f2.Position = f1.Position;

%% SAVE
exportgraphics(f1, "Figures/friction_average_eta="+eta+"_minimum="+B*1e9+"nm.pdf")
exportgraphics(f2, "Figures/friction_directional_eta="+eta+"_minimum="+B*1e9+"nm.pdf")

%% FUNCTIONS
function beta = stokes_drag_sphere(eta, radius)
    beta_r = 6*pi*eta*radius;
    beta = [beta_r, beta_r, beta_r];
end
function beta = stokes_drag_oblate(eta, thickness, diameter)
    e = sqrt(1 - thickness^2/diameter^2);
    a = diameter/2;
    
    beta_ax = 8*pi*eta*a * e^3 / (e*sqrt(1-e^2) - (1 - 2*e^2)*asin(e));
    beta_tr = 16*pi*eta*a * e^3 / ((1 + 2*e^2)*asin(e) - e*sqrt(1-e^2));
    beta = [beta_tr, beta_tr, beta_ax];
end
function beta = stokes_drag_prolate(eta, diameter, length)
    e = sqrt(1 - diameter^2 / length^2);
    a = length/2;
    L = log((1+e)/(1-e));
    beta_ax = 16*pi*eta*a * e^3 / ((1 + e^2)*L - 2*e);
    beta_tr = 32*pi*eta*a * e^3 / (2*e + (3*e^2 - 1)*L);
    beta = [beta_tr, beta_tr, beta_ax];
end
function beta = stokes_drag(eta, diameter, length)
    if diameter == length
        beta = stokes_drag_sphere(eta, diameter/2);
    elseif diameter > length
        beta = stokes_drag_oblate(eta, length, diameter);
    else
        beta = stokes_drag_prolate(eta, diameter, length);
    end
end