function make_plots(psins, plots_params)
	% function make_plots(psins, plots_params)
	% make plots
	% here plots_params = {h, N, s, g0, l, plot_title_str, method_str, plot_ys, plot_zs, plot_energy}.

	h = plots_params{1};
	N = plots_params{2};
	s = plots_params{3};
	g0 = plots_params{4};
	l = plots_params{5};
	plot_title_str = plots_params{6};
	method_str = plots_params{7};
	plot_ys = plots_params{8};
	plot_zs = plots_params{9};
	plot_energy = plots_params{10};

	part_values_str = sprintf('$s=%i$, $h=%.1f$, $N=%i$', s, h, N);
	values_str = sprintf('%s, %s', method_str, part_values_str);

	full_title_str = sprintf('%s\n%s', plot_title_str, values_str);
	t = linspace(0,2*pi);

	% Plot the position coordinates
	figure
	plot(psins(1,:),psins(2,:), 'bo', 'LineWidth', 1.5)
	hold on
	plot(cos(t),sin(t), 'k-', 'LineWidth', 1.5)
	xlabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
	ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
	set(get(gca,'YLabel'),'Rotation',0)
	set(gca,'FontSize',16)
	title(full_title_str,'Interpreter','latex','FontSize',24)
	xlim([-1 1])
	ylim([-1 1])
	axis square
	grid on
	legend('Numerical Solution', 'The Manifold $g(y)=0$', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
	hold off

	if plot_ys
		% Plot y1 vs t, y2 vs t

		titlestr1 = '$y_1$ vs.~$t$';
		titlestr2 = '$y_2$ vs.~$t$';

		figure
		tiledlayout(2,1)
		nexttile
		plot(0:h:N*h, psins(1,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_1$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr1,'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, psins(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$y_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr2,'Interpreter','latex','FontSize',24)
		grid on
		hold off
	end

	if plot_zs
		% plot z1 vs t, z2 vs t
		titlestr3 = '$z_1$ vs.~$t$';
		titlestr4 = '$z_2$ vs.~$t$';

		figure
		tiledlayout(2,1)
		nexttile
		plot(0:h:N*h, psins(3,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$z_1$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr3,'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, psins(4,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$z_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr4,'Interpreter','latex','FontSize',24)
		grid on
		hold off
	end

	if plot_energy
		% make the energy vector
		yns = psins(1:2,:);
		zns = psins(3:4,:);
		ncols = size(psins,2);

		for i=1:ncols
			yn = yns(:,i);
			zn = zns(:,i);
			if i==1
				E = (1/2)*norm(zn)^2 + g0*(l+yn(2));
			else
				Enp = (1/2)*norm(zn)^2 + g0*(l+yn(2));
				E = [E Enp];
			end
		end

		titlestr3 = 'Energy vs.~$t$';
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
end

