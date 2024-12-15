function make_movie(yns,zns,plots_params)
% function make_movie(yns,zns,plots_params)
% make a movie and play it

h = plots_params{1};
N = plots_params{2};

titlestr0 = 'Position Coordinates $(y_1,y_2)$ of the Pendulum';
t = linspace(0,2*pi);

moviein(N);
figure
for j=1:length(yns)
	plot(yns(1,j),yns(2,j), 'bo', 'LineWidth', 1.5)
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
	
	M(:,j) = getframe;
end

% play the movie
movie(M,1,10)


