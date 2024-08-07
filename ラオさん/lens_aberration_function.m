function [phase] = lens_aberration_function(r, f, n)
%LENS_ABERRATION_FUNCTION "thickness" error of aberrated lens.
% use exp(i * k * lens_aberration_function) to modify phase of a field
% Born & Wolf 1964
% https://iopscience.iop.org/article/10.1088/0031-9155/14/2/002/pdf
% L R Evans and C G Morgan 1969 Phys. Med. Biol. 14 205
% r - radial coordinate
% f - lens focal length
% n - lens refractive index

    phase = -0.25 .* (r).^4 .* ...
        (n^2 ./ (8 * f^3 * (n - 1)^2) ...
        - n ./ (8 * (n + 2)) * f ...
        + (n + 1)^2 ./ (2 * n * (n + 2)) * f...
        );
end

