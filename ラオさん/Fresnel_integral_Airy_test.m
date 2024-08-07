%{
test comparison with analytical Airy disk function
%}

clear all

%all parameters in meters
l = 0;
n = 1.5;        % refractive index of lens
f = 30*10^-3;      % lens focal length / m
P = 100;        % laser power / W
a = 1*10^-3;       % aperture radius / m
lam = 694.3*10^-9;   % wavelength / m
k = 2 * pi / lam;     % wavenumber
z = f;          % propagation distance

include_aberration = true;  % select to include aberration term
n_rects = 2000;      % number of trapeziums for manual integration

max_output_radius = 5*10^-6; % m
radius_step_size = max_output_radius / 100; % m

Intensity_inbuilt = [];   % using inbuilt integral function
Intensity_rect = [];      % using manual rectangular integral
Intensity_inbuilt_Seigman = [];   % using inbuilt integral function
Intensity_rect_Seigman = [];      % using manual rectangular integral

Radii = 0:radius_step_size:max_output_radius; % output radii samples

% configure input beam options

% Gaussian profile
input_field_Gauss = @(r) exp(-r.^2 ./ (2 .* a.^2));
input_field_Gauss_xy = @(x, y) exp(-(x.^2 + y.^2) ./ (2 .* a.^2));

% unit intensity radius a
input_field_disk = @(r) 1;
input_field_disk_xy = @(x, y) exp(-((x.^2 + y.^2) ./ (a.^2)).^100);

% phase term to add lens focussing
lens_phase_r = @(r) exp(-1j * k / (2 * f) * r.^2);
lens_phase_xy = @(x, y) exp(-1j * k / (2 * f) * (x.^2 + y.^2));

% aberration term
if include_aberration
    lens_aberration = @(r) exp(1i * k * lens_aberration_function(r, f, n));
else
    lens_aberration = @(r) 1;
end

if include_aberration
    lens_aberration_xy = @(x, y) exp(1i * k * lens_aberration_function(sqrt(x.^2 + y.^2), f, n));
else
    lens_aberration_xy = @(x, y) 1;
end

%% perform different solution routines
disp('Starting Alkelly integral inbuilt')
tic
for R = Radii
    
    % A = integral prefactor
    A = (6.33 * pi * P * a^2 / (lam^2 * f^2));

    % define function with integration parameter r (input radius, r' in eq 9)
    fun = @(r)  (  input_field_disk(r) ...
               .* exp(-1i .* k ./ z ...
                      .* (-Phi(r, f, n, include_aberration) + R .* r ./ 2) ...
                      )...
               .* besselj(0, k .* r .* R ./ z)...
               .* r);
    
    
    E_inbuilt = integral(fun, 0, a);  % evaluated integral
    
    I_inbuilt = A .* (abs(E_inbuilt)).^2;

    Intensity_inbuilt = [Intensity_inbuilt, I_inbuilt];
end
toc

disp("Starting Alkelly integral rect")
tic
for R = Radii

    % A = integral prefactor
    A = (6.33 * pi * P * a^2 / (lam^2 * f^2));

    % define function with integration parameter r (input radius, r' in eq 9)
    fun = @(r)  (  input_field_disk(r) ...
               .* exp(-1i .* k ./ z ...
                      .* (-Phi(r, f, n, include_aberration) + R .* r ./ 2) ...
                      )...
               .* besselj(0, k .* r .* R ./ z)...
               .* r);
    
    E_rect = integral_rect(fun, 0, a, n_rects);
    
    I_rect = A .* (abs(E_rect)).^2;
    
    Intensity_rect = [Intensity_rect, I_rect];
end
toc

disp("Starting Seigman integral inbuilt")
tic
for R = Radii
    % define Seigman function with integration parameter r (input radius, r' in eq 9)
    fun_Seig = @(r)  (  input_field_disk(r) ...
                .* lens_phase_r(r)...
                .* lens_aberration(r)...
               .* exp(-1i .* k ./ z ...
                      .* ( - (r.^2) ./ 2) ...
                      )...
               .* besselj(0, k .* r .* R ./ z)...
               .* r);
    
    E_inbuilt_Seigman = integral(fun_Seig, 0, a);  % evaluated integral
    
    I_inbuilt_Seigman = A .* (abs(E_inbuilt_Seigman)).^2;

    Intensity_inbuilt_Seigman = [Intensity_inbuilt_Seigman, I_inbuilt_Seigman];
end
toc

disp("Starting Seigman integral inbuilt - function")
tic
input_field_func = @(r) (input_field_disk(r) ...
                        .* lens_phase_r(r)...
                        .* lens_aberration(r)...
                        );
E_inbuilt_function = Fresnel_Bessel_integral(0,...
                                            Radii, ...
                                             a, ...
                                             input_field_func, ...
                                             z, ...
                                             lam);
Intensity_inbuilt_function = (abs(E_inbuilt_function)).^2;
toc

disp("Starting Seigman integral rect")
tic
for R = Radii

    % define Seigman function with integration parameter r (input radius, r' in eq 9)
    fun_Seig = @(r)  (  input_field_disk(r) ...
                .* lens_phase_r(r)...
                .* lens_aberration(r)...
               .* exp(-1i .* k ./ z ...
                      .* (- (r.^2) ./ 2) ...
                      )...
               .* besselj(0, k .* r .* R ./ z)...
               .* r);
    E_rect_Seigman = integral_rect(fun_Seig, 0, a, n_rects);
    
    I_rect_Seigman = A .* (abs(E_rect_Seigman)).^2;

    Intensity_rect_Seigman = [Intensity_rect_Seigman, I_rect_Seigman];
end
toc

% Fresnel integral solution


input_field = @(x, y) (sqrt(f^2 * P / (a^2 * pi * (1- exp(-1))))...
                      / f ...
                      * input_field_disk_xy(x, y) ...
                      .* lens_phase_xy(x, y)...
                        .* lens_aberration_xy(x, y)...
                        );

disp("Starting Fresnel integral rect")
tic
E_Fresnel = Fresnel_int(input_field, n_rects, a, lam, z, Radii, zeros(size(Radii)));
%E_Fresnel = 0;   % set to 0 and comment above line to not compute Fresnel_2D integral
toc
Intensity_Fresnel = (abs(E_Fresnel)).^2;

disp("Starting Airy evaluation")
tic
Intensity_Airy = Airy_disk(Radii, 2 * a, f, lam);
toc
%% 
% plot results

Intensity_inbuilt = normalise_intensity(Intensity_inbuilt);
Intensity_rect = normalise_intensity(Intensity_rect);
Intensity_Airy = normalise_intensity(Intensity_Airy);
Intensity_inbuilt_Seigman = normalise_intensity(Intensity_inbuilt_Seigman);
Intensity_rect_Seigman = normalise_intensity(Intensity_rect_Seigman);
Intensity_Fresnel = normalise_intensity(Intensity_Fresnel);
Intensity_inbuilt_function = normalise_intensity(Intensity_inbuilt_function);

close all
hold on

hold on
plot(Radii, Intensity_inbuilt, 'r', 'DisplayName', "inbuilt")
plot(Radii, Intensity_rect, 'b', 'DisplayName', "rect")
plot(Radii, Intensity_Airy, 'g', 'DisplayName', "Airy")
legend
figure
hold on
plot(Radii, Intensity_inbuilt_Seigman, 'r', 'DisplayName', "inbuilt Seigman")
plot(Radii, Intensity_rect_Seigman, 'b', 'DisplayName', "rect Seigman")
plot(Radii, Intensity_Airy, 'g', 'DisplayName', "Airy")
plot(Radii, Intensity_inbuilt_function, 'k', 'DisplayName', "inbuilt function")
legend
figure
hold on
plot(Radii, Intensity_Fresnel, 'r', 'DisplayName', "Fresnel")
plot(Radii, Intensity_Airy, 'g', 'DisplayName', "Airy")
legend

%yline(max(Intensity_rect) / (exp(1)^2), 'DisplayName', "1/e^2 point")

legend

xlim([-inf, inf]);
ylim([0, inf]);

xlabel("Radius / m");
ylabel("Intensity");

function[phase] = Phi(r, f, n, include_aberration)
% phase modifier of aberrated lens
% r - radial coordinate
% f - lens focal length
% n - lens refractive index
% include_aberration - boolean to decide if phase aberrations used

    p = 1 / f;   % lens power
    if include_aberration
            phase = -0.25 .* (r).^4 .* ...
        (n^2 .* p.^3 ./ (8 * f^3 * (n - 1)^2) ...
        - n ./ (8 * (n + 2)) * f ...
        + (n + 1)^2 ./ (2 * n * (n + 2)) * f...
        );
    else
        phase = 0;
    end
end

function[integral] = integral_rect(func, x_min, x_max, n_rects)
% perform trapzeium rule integration
% func(x) - functional form of function to integrate
% x_min - lower bound of integral
% x_max - upper bound of integral
% n_rects - number of trapeziums to use

dx = (x_max - x_min) / n_rects;

x = linspace(x_min, x_max, n_rects);

integral = dx * sum(  func( x(2:(n_rects-1)) )  ) ...
         + dx / 2 * (  func(x(1)) + func(x(n_rects))  );
end

function[integral] = integral2D_rect(f, x, y)
% array inputs
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

function[field] = Fresnel_int(source_field, source_step_number, a, wavelength, z, output_x, output_y)
% compute the Fresnel 2D integral for a given source field and distance
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
        field(i, j) = integral2D_rect(U0 .* exp( 1i * k / (2 * z) ...
                                                * ( (x(i, j) - x0).^2 ...
                                                    + (y(i, j) - y0).^2 ...
                                                  )...
                                               )...
                                      , x0, y0);
    end
end

field = field * exp(1i * k * z) / (1i * wavelength * z);
end

function[intensity] = Airy_disk(radius, diameter, focal_length, wavelength)
% the Airy disk function - a unit circle field brought to focus by a lens
% radius = radius coordinates to compute
% diameter = diameter of unit disk input to lens
% focal_length = focal length of lens
% wavelength = wavelength

D = diameter;
f = focal_length;
w = wavelength;

NA = D / (2 * f);   % numerical aperture

intensity = abs(...
                besselj(1, 2 * pi * radius * NA / w) ...
                ./ (2 * pi * radius * NA / w)...
                ).^2;

intensity(radius == 0) = 1/4;  % fix r=0 coordinate

end

function[intensity] = normalise_intensity(intensity)
% normalise intensity to peak=1 of distribution

intensity = intensity ./ max(max(intensity));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%