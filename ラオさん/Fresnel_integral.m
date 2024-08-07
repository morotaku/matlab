function[field] = Fresnel_integral(source_field, source_step_number, a, wavelength, z, output_x, output_y)
% FRESNEL_INTEGRAL compute the Fresnel 2D integral for a given source field and distance
% source_field(x, y) = functional form of source input
% source_step_number = number of samples of source field to take
% a = aperture (integration) radius
% wavelength
% z = propagation distance
% output_radius = output field grid radius
% output_x = output x coordinates (fleshed out)
% output_y = output y coordinates (fleshed out)

k = 2 * pi / wavelength;

x0 = linspace(-a, a, source_step_number);
[x0, y0] = meshgrid(x0, x0);

U0 = source_field(x0, y0);

x = output_x;
y = output_y;

output_size = size(x);
field = zeros(output_size);

for i = 1:output_size(1)
    for j = 1:output_size(2)
        field(i, j) = integral2D_trapezium(U0 .* exp( 1i * k / (2 * z) ...
                                                * ( (x(i, j) - x0).^2 ...
                                                    + (y(i, j) - y0).^2 ...
                                                  )...
                                               )...
                                      , x0, y0);
    end
end

field = field * exp(1i * k * z) / (1i * wavelength * z);
end

function[integral] = integral2D_trapezium(f, x, y)
% perform 2D trapezium integration on function
% f = array of function values
% x = x coordinate grid
% y = y coordinate grid

temp = size(x);
array_size = temp(1);

dx = (max(max(x)) - min(min(x))) / array_size;

integral = dx^2 * sum(sum(  f(2:(array_size-1), 2:(array_size-1) )  )) ...
         + dx^2 / 4 * (  f(1, 1) + f(1, array_size) + f(array_size, 1) + f(array_size, array_size)  )...
         + dx^2 / 2 * (...
                          sum(sum(f(1, 2:(array_size-1)))) ...
                        + sum(sum(f((array_size-1), 2:(array_size-1)))) ...
                        + sum(sum(f(2:(array_size-1), 1))) ...
                        + sum(sum(f(2:(array_size-1), (array_size-1))))  ...
                        );
end