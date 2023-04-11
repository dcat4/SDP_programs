function [cv_ind] = get_kfold_cv_indices(data, k)
% creates k x size(data,1)/k matrix of random indices that can be applied
% to subset data for k-fold cross-validation

rand_ns = randperm(size(data,1)); % all the random numbers you need
cv_ind = zeros(k , ceil(size(data,1) / k));
n_leftovers = ( (size(data,1) / k) - ( floor(size(data,1) / k) )) * k;
counter_start = n_leftovers + 1;
counter_end = n_leftovers + floor(size(data,1)/k);
for i = 1:k
    cv_ind(i,1:floor( size(data,1) / k)) = rand_ns(counter_start:counter_end);
    counter_start = counter_start + floor(size(data,1)/k);
    counter_end = counter_end + floor(size(data,1)/k);
end
% add the leftovers to the CV_indices array, and put NaN's for the sets
% where there's not enough leftovers to go around
leftovers = rand_ns(1:n_leftovers);
na_array = ones(1, k - length(leftovers)) .* NaN;
leftovers = [leftovers , na_array];
cv_ind(:,length(cv_ind)) = leftovers;

