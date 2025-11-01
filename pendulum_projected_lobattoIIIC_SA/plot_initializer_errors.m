function [alpha beta] = plot_initializer_errors(stepsizes, errs, plot_errs_params)
% function [alpha beta] = plot_initializer_errors(stepsizes, errs, plot_errs_params)
% Finds the least squares line log10(errs) = alpha + beta*log10(H).
%
% I.e., fits the least squares regression line to the data points (x_i,y_i),
% where x_i = log10(errs(i)) and y_i = log10(stepsizes(i)) and returns [alpha, beta],
% where alpha is the intercept and beta is the slope.
%
% Also plots the least-squares line in log-log space with title given
% by title_str.
%
% Also prints the equation of the line using print_str.
%
% Inputs:
%	stepsizes: vector of stepsize values (just h, not rh)
%	errs: vector of errors in starting approximations
% 	plot_errs_params = {title_str, print_str, r}: parameters 
%	title_str: string used in the title
%	print_str: string used in the print statment
%
% Outputs:
%	alpha: intercept of the least squares regression line
%	beta: slope of the least squares regression line

%% unpack parameters
title_str = plot_errs_params{1};
print_str = plot_errs_params{2};

%% find the least squares line in log-log space
x = log10(stepsizes);
y = log10(errs);

% make the design matrix
A = ones(length(x),2);
A(:,2) = x;

% Solve for the coefficients
coeff = (A' * A) \ (A' * y);
alpha = coeff(1);
beta = coeff(2);

x_reg = linspace(min(log10(stepsizes)),max(log10(stepsizes)));
y_reg = alpha + beta * x_reg;

%% print the line
fprintf('\nLeast squares line for ')
fprintf(print_str)
fprintf(' : log10(errs) = %0.2f + %0.2f * log10(stepsizes)\n', alpha, beta)

eq_str = sprintf('Least Squares Regression Line: \n');
eq_str = append(eq_str,'$\log_{10}($errs$) = $');
eq_str = append(eq_str, sprintf('%0.2f + %0.2f',alpha,beta));
eq_str = append(eq_str, '$\times \log_{10}(h)$');

%% plot the errors and the least-squares line
figure
plot(x, y,'ro--','LineWidth',2.0)
hold on
plot(x_reg, y_reg, 'b-', 'LineWidth', 2.0)
xlabel('$\log_{10}(h)$','Interpreter','latex','FontSize',22)
ylabel('$\log_{10}($Error$)$','Interpreter','latex','FontSize',22)
set(gca,'FontSize',20)
title(title_str,'Interpreter','latex','FontSize',28)
legend('Error',eq_str,'Interpreter','latex','FontSize',22,'Location','Northwest')
ylim_current  = ylim;
new_ticks = linspace(ylim_current(1), ylim_current(2), 3);
yticks(new_ticks);
hold off


