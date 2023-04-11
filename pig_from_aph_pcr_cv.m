function [coefficients, intercepts, summary_gofs, all_gofs, pig_predicted] = pig_from_aph_pcr_cv(aph, pig,...
    n_cv, cv_frac, max_components, pred_select, varargin)
% builds and cross-validates principal component regression models to
% predict pigment concentrations from aph (or any optical spectra)
% aph is an observation x variable matrix
% pig is a vector (only 1 pigment can be modeled at a time)
% n_cv is the number of cross-validations to perform
% cv_frac is the fraction of training data for each cv
% max_components is the max number of pcs that can be included in pcr
% models
% pred_select is either 'stepwise' or 'kfold descent' to specify the method
% for selecting/optimizing regression models
% varargin is either the 'Criterion' for stepwise regression model 
% selection using stepwiselm (see matlab documentation) or the
% 'mdl_pick_metric' for my custom kfold descent algorithm

    % set the random number generator so you're random cross-validating is
    % reproducible...
    rng(1); 
    
    % the CV loop:
    for i = 1:n_cv
        % Create broad training data (75%) and validation data (25%)
        % do split_data_2way_random on pigments and spectra
        [aph_train, aph_val, idx] = split_data_2way_random(aph, cv_frac);
        pigs_train = pig(idx); 
        pigs_val = pig; pigs_val(idx) = [];

        % optimize model coefficients on training data:
        if strcmp(pred_select, 'stepwise')

            [betas, alpha] = find_best_model(aph_train, pigs_train, ...
                pred_select, varargin{1}, max_components);

        elseif strcmp(pred_select, 'kfold descent')

            [betas, alpha] = find_best_model(aph_train, pigs_train, ...
                pred_select, 5, varargin{1}, max_components);
            % k is hard-coded at 5

        end

        % Store mean/std of each set of k-fold CV betas (the model coefficients for the ith run of the n_permutations):
        mean_betas(:,i) = mean(betas , 2);
        mean_alphas(i) = mean(alpha);
        std_betas(:,i) = std(betas,0,2);
        std_alphas(i) = std(alpha);

        % Validate on the data you set aside previously:
        modeled_pigs = (aph_val * mean_betas(:,i)) + mean_alphas(i);

        % constrain modeled pigments:
        modeled_pigs(modeled_pigs < 0) = 0;
        
        % this cv's validation statistics:
        lin_mdl_4_validation = fitlm(modeled_pigs, pigs_val);
        R2s_final(i) = lin_mdl_4_validation.Rsquared.Ordinary;
        RMSEs_final(i) = lin_mdl_4_validation.RMSE;

        pigs_val(pigs_val == 0) = 1e-4; % add a really small # to zeros to make % error calculations reasonable (pft_index-independent)
        % calc all the validation statistics
        pct_bias(i) = mean(((modeled_pigs - pigs_val)./pigs_val)*100);
        pct_errors(i,:) = abs(((modeled_pigs - pigs_val)./pigs_val)*100);
        med_pct_error(i) = median(pct_errors(i,:));
        avg_pct_error(i) = mean(pct_errors(i,:));
        sort_pct_errors = sort(pct_errors(i,:));
        CI_pct_error(i) = sort_pct_errors(ceil((0.95 * size(pct_errors(i,:),2))));
        std_pct_error(i) = std(pct_errors(i,:));

        % save predicted pigments:
        pig_predicted(i, :) = modeled_pigs;
        
        % display permutation number here...
        disp(['on CV number ', num2str(i)]);

    end
    
    % compile outputs and save them appropriately.
    coefficients = mean_betas;
    intercepts = mean_alphas;
    
    % get summary gof stats and pop them in a table:
    summary_gofs = [mean(R2s_final),std(R2s_final),mean(RMSEs_final),std(RMSEs_final),mean(avg_pct_error),std(avg_pct_error),...
        mean(med_pct_error),std(med_pct_error),mean(pct_bias),std(pct_bias)];
    summary_gofs = array2table(summary_gofs);
    summary_gofs.Properties.VariableNames = {'Mean_R2','SD_R2','Mean_RMSE','SD_RMSE',...
        'Mean_mean_pct_error','SD_mean_pct_error','Mean_median_pct_error','SD_median_pct_error',...
        'Mean_pct_bias','SD_pct_bias'};
    
    % get all gof stats from n_permutations CV's and pop them into a struct
    all_gofs = table;
    all_gofs.R2s = R2s_final;
    all_gofs.RMSEs = RMSEs_final;
    all_gofs.mean_pct_error = avg_pct_error;
    all_gofs.median_pct_error = med_pct_error;
    all_gofs.pct_bias = pct_bias; 
    
end