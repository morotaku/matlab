function [field] = Fresnel_Bessel_integral(azimuthal_symmetry_order,...
                                           output_radii, ...
                                          input_aperture_radius, ...
                                          input_field, ...
                                          propagation_distance, ...
                                          wavelength)

%FRESNEL_BESSEL_INTEGRAL compute Fresnel diffraction with cylindrical symmetry
% azimuthal_symmetry_order = (m) integer, exp((+-) m theta) symmetry of
%                            input field
% output_radii = array of output radius values to evaluate (meters)
% input_aperture_radius = radius of input field to integrate to (meters)
% input_field = the input field in functional form
% propagation_distance = distance from input to output planes (meters)
% wavelength = wavelength of field (meters)
%
% Equation source: Eq.(11)
% 10.1364/AO.32.005893
% Hu, Y., Mo, G., Ma, Z., Fu, S., Zhu, S., Yin, H., Li, Z., Chen, Z.
% (2021). Vector vortex state preservation in Fresnel cylindrical
% diffraction. Optics Letters, 46(6), 1313.

k = 2 * pi / wavelength;
z = propagation_distance;

field = [];

for R = output_radii

    integrand = @(r)  (   input_field(r)...
                       .* exp(1i .* k ./ z .* ( (r.^2) ./ 2))...
                       .* besselj(azimuthal_symmetry_order, k .* r .* R ./ z)...
                       .* r);

    field_sample = integral(integrand, 0, input_aperture_radius);  % evaluate integral;

    field = [field, field_sample];
end

% add integral prefactors
field = field .* (k * ((-1j) ^ (1 + azimuthal_symmetry_order)) * exp(1j * k * z) / z...
                                .* exp(1j * k * output_radii.^2 / (2 * z))...
                                );
end
