%{
plot x-z slice of beam under spherical aberration
%}

%% configure simulation parameters
clear all
close all
clc
% input beam parameters
wav = 640*10^-9;   % wavelength / m
input_beam_radius = 2*10^-3; %  / m
l = 0;  % exp(i m theta) symmetry of input

% lens parameters
f = 25.4*10^-3;      % lens focal length / m
n = 1.5;        % refractive index of lens
include_aberration = true;  % select to include aberration term

% simulation parameters
integration_radius = 3 * sqrt(1 + abs(l)) * input_beam_radius;       % / m
max_output_radius = 100*10^-6; % m
output_radius_sample_number = 200;

distance_half_range = 2*0.5*1e-3; % /m
distance_start = f - distance_half_range;
distance_end = f + distance_half_range/10;   %  / m
distance_sample_number = 200;

%% set up simulation

k = 2 * pi / wav;     % wavenumber

Radii = linspace(0, max_output_radius, output_radius_sample_number); % output radii samples
Radii_step = Radii(2) - Radii(1);

z = linspace(distance_start, distance_end, distance_sample_number); % output distance sample
z_step = z(2) - z(1);

[Radii, z] = meshgrid(Radii, z);

% Gaussian profile
input_field_Gauss = @(r) r .^ abs(l) .* exp(-r.^2 ./ (2 .* input_beam_radius.^2));
%input_field_Gauss = @(r) r .* exp(-r.^2 ./ (2 .* input_beam_radius.^2));
%input_field_Gauss = @(r) exp(-(r.^2 ./ (2 .* input_beam_radius.^2)).^10);

% phase term to add lens focussing
lens_phase_r = @(r) exp(-1j * k / (2 * f) * r.^2);

% aberration term
lens_aberration = @(r) exp(1i * k * lens_aberration_function(r, f, n));

% build whole function
if include_aberration
    input_field = @(r) (input_field_Gauss(r) .* lens_phase_r(r) .* lens_aberration(r));
else
    input_field = @(r) (input_field_Gauss(r) .* lens_phase_r(r));
end

output_field = zeros(size(Radii));

for i = 1:distance_sample_number
    output_field(i, :) = Fresnel_Bessel_integral(l,...
                                                 Radii(i, :), ...
                                                integration_radius, ...
                                                input_field, ...
                                                z(i, 1), ...
                                                wav);
end

output_intensity = (abs(output_field)).^2;

%% plot results

%% check power is conserved

r_samples = linspace(0, integration_radius, 100);
r_samples_step = r_samples(2) - r_samples(1);

input_power = sum(2 * pi * r_samples .* (abs(input_field(r_samples)).^2) * r_samples_step);

power_check = [];
for i = 1:distance_sample_number
        power_check = [power_check, sum(2 * pi * Radii(i, :) .* output_intensity(i, :)) * Radii_step];
end
power_check = power_check / input_power;

figure(1)
clf

plot(z(:, 1), power_check)
title('check simulation power conservation')
xlabel('distance / m')
ylabel('ratio beam power to input power')
xlim([distance_start, distance_end])

%% plot intensity 2D

% duplicate +ve radii to -ve
intensity_plot = [output_intensity(:, end:-1:2) , output_intensity];
Radii_plot = [-Radii(:, end:-1:2) , Radii];
z_plot = [z(:, end:-1:2) , z];

figure(2)
clf

surf(z_plot, Radii_plot, intensity_plot, 'EdgeColor', 'None')
title('intensity')
xlabel('distance from lens / m')
ylabel('radius / m')
view(2)
xlim([distance_start, distance_end])
ylim([-max_output_radius, max_output_radius])

% plot field magnitude 2D
figure(3)
clf

surf(z_plot, Radii_plot, sqrt(intensity_plot), 'EdgeColor', 'None')
title('field magnitude')
xlabel('distance from lens / m')
ylabel('radius / m')
view(2)
xlim([distance_start, distance_end])
ylim([-max_output_radius, max_output_radius])

% plot field intensity normalised to peak
figure(4)
clf

intensity_normalised = intensity_plot;
for i = 1:distance_sample_number
    intensity_normalised(i, :) = intensity_normalised(i, :) / max(max(intensity_normalised(i, :)));
end

surf(z_plot, Radii_plot, intensity_normalised, 'EdgeColor', 'None')
title('intensity normalised')
xlabel('distance from lens / m')
ylabel('radius / m')
view(2)
xlim([distance_start, distance_end])
ylim([-max_output_radius, max_output_radius])

%% plot intensity cross sections
z_samples = [0.0298, 0.0298 - 200e-6, 0.0298 + 200e-6];    % select which z distances to plot

figure(5)
clf
hold on

for z_sample = z_samples
    [value, index] = min(abs(z_plot(:, 1) - z_sample));
    plot(Radii_plot(index, :), intensity_plot(index, :), 'DisplayName', 'z = ' + num2str(z_sample))
end

xlabel('radius / m')
ylabel('intensity / arb.')

hold off
legend
