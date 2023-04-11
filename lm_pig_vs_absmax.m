function [R2, amax_wl] = lm_pig_vs_absmax(pigmat, piglabels, amat, wl)
% fits linear regression model probe the relationship between sums of
% pigment concentrations and their known absorption maxima
% pigmat is a matrix of pigment concentrations
% amat is a matrix of (smoothed) absorption spectra
% R2 is a vector of R2

% from the literature:
pigXamax = {["TChla"], 437; ...
    ["MVChlb" , "Chlc1c2", "Chlc3", "Zea", "Diadino", "Diato", "Allo", "ABCar"], 465;...
    ["Zea", "Diadino", "Diato", "Allo", "ABCar"], 492; ...
    ["Fuco", "Perid"], 542; ...
    ["Chlc1c2", "Chlc3"], 590};

% get 2nd derivative of input spectra:
[~, d2, dwl] = calc_spec_derivative(amat, wl);

% do regression:
for i = 1:size(pigXamax,1)

    % get sum of pigments:
    pigi = pigXamax{i,1};
    pigsum = sum(pigmat(:, ismember(piglabels, pigi)), 2);

    % get max absorption:
    amax_wl = pigXamax{i,2};
    amax = d2(:, dwl == amax_wl);

    % fit regression model:
    lm = fitlm(amax, pigsum);
    R2(i) = lm.Rsquared.Adjusted;
end

