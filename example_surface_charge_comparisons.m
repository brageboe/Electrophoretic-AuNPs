clear variables 
A = linspace(2, 200, 1e3).*1e-9;  % length of one gnp axis (diameter or length)
B = min(A);     % length of other gnp axis (length or diameter)
n = 5;
e0 = 1.602176634e-19;

% SPHERES
for i=1:length(A)
    q1(i) = net_charge_Gupta(A(i),A(i));
    q2(i) = net_charge_Rostalski(A(i),A(i));
end
f = figure();
subplot(1,3,1)
semilogy(A/1e-9, q1/e0)
hold on
semilogy(A/1e-9, q2/e0)
title("Spheres")
xlabel("Diameter (nm)")
ylabel("Surface charge ($e_0$)")
ylim([1e0 1e6])
grid on
hold off

% PROLATES
for i=1:length(A)
    q1(i) = net_charge_Gupta(B,A(i));
    q2(i) = net_charge_Rostalski(B,A(i));
end
subplot(1,3,2)
semilogy(A/1e-9, q1/e0)
hold on
semilogy(A/1e-9, q2/e0)
title("Prolates")
xlabel("Length (nm)")
ylim([1e0 1e6])
grid on
hold off

% OBLATES
for i=1:length(A)
    q1(i) = net_charge_Gupta(A(i),B);
    q2(i) = net_charge_Rostalski(A(i),B);
end
subplot(1,3,3)
semilogy(A/1e-9, q1/e0)
hold on
semilogy(A/1e-9, q2/e0)
title("Oblates")
xlabel("Diameter (nm)")
ylim([1e0 1e6])
grid on
hold off

axesHandles = get(gcf, 'children');
set(axesHandles, 'LineWidth', 1);
set(axesHandles, 'FontSize', 14);
set(axesHandles(1), 'YTickLabel', []);
set(axesHandles(2), 'YTickLabel', []);
legend("Gupta (n="+n+")", "Rostalski", 'Location', 'SouthEast', 'FontSize', 10)

% exportgraphics(f,"Figures/surfacecharge_comparisons_minimum="+B*1e9+"nm.png")