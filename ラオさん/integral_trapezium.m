function [integral] = integral_trapezium(func, x_min, x_max, n_rects)
% INTEGRAL_TRAPEZIUM perform trapzeium rule integration
% func(x) - functional form of function to integrate
% x_min - lower bound of integral
% x_max - upper bound of integral
% n_rects - number of trapeziums to use

dx = (x_max - x_min) / n_rects;

x = linspace(x_min, x_max, n_rects);

integral = dx * sum(  func( x(2:(n_rects-1)) )  ) ...
         + dx / 2 * (  func(x(1)) + func(x(n_rects))  );
end

