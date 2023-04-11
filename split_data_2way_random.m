function [split1, split2, idx] = split_data_2way_random(data, frac1)
% splits an observation X variable matrix into 2 random subsets with split1
% containing frac1 of the total observations
% idx provides entries of data that are included in split 1

idx = randperm(size(data,1) , floor(size(data,1) * frac1))'; 
split1 = data(idx , :); 
split2 = data; 
split2(idx , :) = []; 
