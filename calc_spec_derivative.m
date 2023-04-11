function [d1, d2, dwl] = calc_spec_derivative(spectra, wl)
% calculates 1st and 2nd derivatives of spectral data using a 2nd order
% centered finite difference approximation.
% assumes input spectra are sampled uniformly
% input spectra is an observation x spectra matrix. wl is a vector of
% length size(spectra,2)
% outputs d1 and d2 are 1st and 2nd derivative spectra along wavelengths
% dwl

band_sep = wl(2) - wl(1); % band separation

d1 = (spectra(: , 3:end) - spectra(: , 1:end-2)) ./ (2 * band_sep); % 1st derivative

d2 = (spectra(: , 3:end) + spectra(: , 1:end-2) - (2 .* spectra(: , 2:end-1)) ) ./ (band_sep ^ 2); % 2nd derivative

dwl = wl(2:end-1); % same for first and second derivative

% remove wl's with nan:
rmme = all(isnan(d2), 1);
d1 = d1(:, ~rmme);
d2 = d2(:, ~rmme);
dwl = dwl(~rmme);
