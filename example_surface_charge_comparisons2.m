% Same as 'example_surface_charge_comparisons.m' but without Gupta charge
clear variables 
A = linspace(3, 200, 1e3).*1e-9;  % length of one gnp axis (diameter or length)
B = min(A);     % length of other gnp axis (length or diameter)
e0 = 1.602176634e-19;

% SPHERES
for i=1:length(A)
    q(i) = net_charge_Rostalski(A(i),A(i));
end
f = figure();
subplot(1,3,1)
semilogy(A/1e-9, q/e0)
hold on
title("Spheres")
xlabel("Diameter [nm]")
ylabel("Surface charge $q$ [$e_0$]")
ylim([1e0 1e4])
grid on
hold off

% PROLATES
for i=1:length(A)
    q(i) = net_charge_Rostalski(B,A(i));
end
subplot(1,3,2)
semilogy(A/1e-9, q/e0)
hold on
title("Prolates")
xlabel("Length [nm]")
ylim([1e0 1e4])
grid on
hold off

% OBLATES
for i=1:length(A)
    q(i) = net_charge_Rostalski(A(i),B);
end
subplot(1,3,3)
semilogy(A/1e-9, q/e0)
hold on
title("Oblates")
xlabel("Diameter [nm]")
ylim([1e0 1e4])
grid on
hold off

% set_standard_plot_style();

axesHandles = get(gcf, 'children');
set(axesHandles, 'LineWidth', 1);
set(axesHandles, 'FontSize', 12);
set(axesHandles(1), 'YTickLabel', []);
set(axesHandles(2), 'YTickLabel', []);

f.Position = [680 458 560 280]; % [x0 y0 width height]
% exportgraphics(f,"Figures/surfacecharge_comparisons2_minimum="+B*1e9+"nm.pdf")