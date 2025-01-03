function make_plots(yns, plots_params)
	% function make_plots(yns, plots_params)
	% make plots (two body system)
	% here plots_params = {h, N, s, e, alpha0, plot_title_str, method_str, plot_qs, plot_ps, plot_energy}.

	h = plots_params{1};
	N = plots_params{2};
	s = plots_params{3};
	e = plots_params{4};
	alpha0 = plots_params{5};
	plot_title_str = plots_params{6};
	method_str = plots_params{7};
	plot_qs = plots_params{8};
	plot_ps = plots_params{9};
	plot_energy = plots_params{10};

	values_str = sprintf('%s, $s=%i$, $h=%.1f$, $N=%i$', method_str, s, h, N);

	full_title_str = sprintf('%s\n%s', plot_title_str, values_str);

	% Plot the position coordinates
	figure
	plot(yns(1,:),yns(2,:), 'b-', 'LineWidth', 1.5)
	hold on

	% plot the analytical solution
	theta = linspace(0,2*pi);
	radius = alpha0 ./ (1 + e*cos(theta));
	x1_true = radius .* cos(theta);
	x2_true = radius .* sin(theta);
	plot(x1_true,x2_true,'r--','LineWidth',2)

	xlabel('$q_1$', 'Interpreter', 'latex', 'FontSize', 16)
	ylabel('$q_2$', 'Interpreter', 'latex', 'FontSize', 16)
	set(get(gca,'YLabel'),'Rotation',0)
	set(gca,'FontSize',16)
	title(full_title_str,'Interpreter','latex','FontSize',24)
	xlim([-2.2 1.2])
	ylim([-1.5 2])
	xticks([-2 -1 0 1])
	yticks([-1 0 1])
	axis square
	grid on
	legend('Numerical Solution', 'True Solution', 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'NorthEast')
	hold off

	if plot_qs
		% Plot q1 vs t, q2 vs t

		titlestr1 = '$q_1$ vs.~$t$';
		titlestr2 = '$q_2$ vs.~$t$';

		figure
		tiledlayout(2,1)
		nexttile
		plot(0:h:N*h, yns(1,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$q_1$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr1,'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, yns(2,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$q_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr2,'Interpreter','latex','FontSize',24)
		grid on
		hold off
	end

	if plot_ps
		% plot p1 vs t, p2 vs t
		titlestr3 = '$p_1$ vs.~$t$';
		titlestr4 = '$p_2$ vs.~$t$';

		figure
		tiledlayout(2,1)
		nexttile
		plot(0:h:N*h, yns(3,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$p_1$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr3,'Interpreter','latex','FontSize',24)
		grid on
		hold off

		nexttile
		plot(0:h:N*h, yns(4,:), 'b-', 'LineWidth', 1.5)
		hold on
		xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16)
		ylabel('$p_2$', 'Interpreter', 'latex', 'FontSize', 16)
		set(get(gca,'YLabel'),'Rotation',0)
		set(gca,'FontSize',16)
		title(titlestr4,'Interpreter','latex','FontSize',24)
		grid on
		hold off
	end

	if plot_energy
		% Make the energy vector.
		% Note that here energy E = H(q,p).
		qns = yns(1:2,:);
		pns = yns(3:4,:);
		ncols = size(yns,2);

		for i=1:ncols
			qn = qns(:,i);
			pn = pns(:,i);
			if i==1
				E = (1/2)*norm(pn)^2 - 1/norm(qn);
			else
				Enp = (1/2)*norm(pn)^2 - 1/norm(qn);
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

