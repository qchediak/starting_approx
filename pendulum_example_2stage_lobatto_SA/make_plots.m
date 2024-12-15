function make_plots(yns,zns,energy_vec,plots_params)
% function make_plots(yns,zns,energy_vec,plots_params)
% make plots

r = plots_params{1};
h = plots_params{2};

t_vals = [0 h (r+1)*h];

E = energy_vec;

N = 2; % number of steps

% Plot the position coordinates
titlestr0 = sprintf('Position Coordinates $(y_1,y_2)$ of the Pendulum, $r$=%0.1f',r);
t = linspace(0,2*pi);

figure
plot(yns(1,:),yns(2,:), 'bo', 'LineWidth', 1.5)
hold on
plot(cos(t),sin(t), 'k-', 'LineWidth', 1.5)
xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
set(get(gca,'YLabel'),'Rotation',0)
set(gca,'FontSize',16)
title(titlestr0,'Interpreter','latex','FontSize',24)
xlim([-1 1])
ylim([-1 1])
axis square
grid on
legend('Numerical Solution', 'The Manifold $g(y)=0$', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'NorthEast')
hold off

% Plot y1 vs t, y2 vs t, and energy vs t

titlestr1 = sprintf('$y_1$ vs.~$t$, $r$=%0.1f',r);
titlestr2 = sprintf('$y_2$ vs.~$t$, $r$=%0.1f',r);
titlestr3 = sprintf('Energy vs.~$t$, $r$=%0.1f',r);

figure
tiledlayout(3,1)
nexttile
plot(t_vals, yns(1,:), 'b-', 'LineWidth', 1.5)
hold on
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
set(get(gca,'YLabel'),'Rotation',0)
set(gca,'FontSize',16)
title(titlestr1,'Interpreter','latex','FontSize',24)
grid on
hold off

nexttile
plot(t_vals, yns(2,:), 'b-', 'LineWidth', 1.5)
hold on
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
set(get(gca,'YLabel'),'Rotation',0)
set(gca,'FontSize',16)
title(titlestr2,'Interpreter','latex','FontSize',24)
grid on
hold off

nexttile
plot(t_vals, E, 'b-', 'LineWidth', 1.5)
hold on
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Energy', 'Interpreter', 'latex', 'FontSize', 16)
set(get(gca,'YLabel'),'Rotation',0)
set(gca,'FontSize',16)
title(titlestr3,'Interpreter','latex','FontSize',24)
grid on
hold off

% Plot energy vs. t in a separate figure
figure
plot(t_vals, E, 'b-', 'LineWidth', 1.5)
hold on
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Energy', 'Interpreter', 'latex', 'FontSize', 16)
set(get(gca,'YLabel'),'Rotation',0)
set(gca,'FontSize',16)
title(titlestr3,'Interpreter','latex','FontSize',24)
grid on
hold off
end

