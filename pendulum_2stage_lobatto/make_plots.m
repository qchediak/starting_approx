function make_plots(yns,zns,energy_vec,plots_params)
	% function make_plots(yns,zns,energy_vec,plots_params)
	% make plots

	h = plots_params{1};
	N = plots_params{2};
	s = plots_params{3};

	E = energy_vec;

	line2 = sprintf('Lobatto IIIA-IIIB, $s=%i$, $h=%0.1f$, $N=%i$', s, h, N);

	% Plot the position coordinates
	titlestr0 = 'Position Coordinates $(y_1,y_2)$ of the Pendulum';
	g_title_str = sprintf('%s\n%s', titlestr0, line2);
	t = linspace(0,2*pi);

	figure
	plot(yns(1,:),yns(2,:), 'bo', 'LineWidth', 1.5)
	hold on
	plot(cos(t),sin(t), 'k-', 'LineWidth', 1.5)
	xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
	ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
	set(get(gca,'YLabel'),'Rotation',0)
	set(gca,'FontSize',16)
	title(g_title_str,'Interpreter','latex','FontSize',24)
	xlim([-1 1])
	ylim([-1 1])
	axis square
	grid on
	legend('Numerical Solution', 'The Manifold $g(y)=0$', 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'NorthEast')
	hold off

	% plot 2y^Tz vs. time
	for i=1:length(yns)
		if i==1
			values = 2*dot(yns(:,i), zns(:,i));
		else
			values = [values 2*dot(yns(:,i), zns(:,i))];
		end
	end

	line1 = 'Adherence of the Solution to $(g_y f)(y,z)=0$';

	gyf_title_str = sprintf('%s\n%s', line1, line2);

	figure
	plot(0:h:N*h, values, 'b-', 'LineWidth', 1.5)
	hold on
	%plot(0:h:N*h, zeros(1,N+1), 'r-', 'LineWidth', 1.5)
	xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
	ylabel('$2y(t)^Tz(t)$', 'Interpreter', 'latex', 'FontSize', 16)
	set(get(gca,'YLabel'),'Rotation',0)
	set(gca,'FontSize',16)
	title(gyf_title_str,'Interpreter','latex','FontSize',24)
	ylim([-0.2 0.2])
	axis square
	grid on
	legend('Numerical Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
	hold off

	% Plot y1 vs t, y2 vs t, and energy vs t

	titlestr1 = '$y_1$ vs.~$t$';
	titlestr2 = '$y_2$ vs.~$t$';
	titlestr3 = 'Energy vs.~$t$';

	figure
	tiledlayout(3,1)
	nexttile
	plot(0:h:N*h, yns(1,:), 'b-', 'LineWidth', 1.5)
	hold on
	xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
	ylabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
	set(get(gca,'YLabel'),'Rotation',0)
	set(gca,'FontSize',16)
	title(titlestr1,'Interpreter','latex','FontSize',24)
	grid on
	hold off

	nexttile
	plot(0:h:N*h, yns(2,:), 'b-', 'LineWidth', 1.5)
	hold on
	xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
	ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
	set(get(gca,'YLabel'),'Rotation',0)
	set(gca,'FontSize',16)
	title(titlestr2,'Interpreter','latex','FontSize',24)
	grid on
	hold off

	nexttile
	plot(0:h:N*h, E, 'b-', 'LineWidth', 1.5)
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
	plot(0:h:N*h, E, 'b-', 'LineWidth', 1.5)
	hold on
	xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
	ylabel('Energy', 'Interpreter', 'latex', 'FontSize', 16)
	set(get(gca,'YLabel'),'Rotation',0)
	set(gca,'FontSize',16)
	title(titlestr3,'Interpreter','latex','FontSize',24)
	grid on
	hold off

end

