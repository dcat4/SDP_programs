function [betas, alpha] = find_best_model(x, y, how, varargin)
% finds the best (linear) model for predicting y from x and returns the
% coefficients (betas) and intercept (alpha)
% how is 'stepwise' or 'kfold descent'
% stepwise can be used with varargin to specify criteria (see matlabs
% stepwiselm function 'Criterion') 
% kfold descent does my custom approach. varargin includes k (a scalar value),
% mdl_pick_metric (a string to say how to pick the best model), and 
% max_components (a scalar value).

if strcmp(how, 'stepwise')
   
    max_components = varargin{2};
    x_z = zscore(x);
    [EOFs,AFs,~,~,~] = pca(x_z,'Centered',false,'NumComponents',max_components,'Rows','complete');

    lm = stepwiselm(AFs, y, 'Upper','linear', 'Criterion', varargin{1}, 'Verbose', 0);
    betas = table2array(lm.Coefficients(:,1));
    alpha = betas(1); % intercept
    betas(1) = [];
    idx = strrep(lm.PredictorNames, 'x', '');
    idx = cellfun(@str2num, idx);

    % in z-scored derivative space:
    betas = EOFs(: , idx) * betas;

    % un-zscore coefficients:
    alpha = alpha - sum(betas .* (mean(x)' ./ std(x)'));
    betas = betas ./ std(x)';

elseif strcmp(how, 'kfold descent')

    % define k, mdl_pick_metric, max_components from varargin:
    k = varargin{1}; 
    mdl_pick_metric = varargin{2};
    max_components = varargin{3};

    % get CV_indices for training:
    CV_indices = get_kfold_cv_indices(x, k);

    % k-fold cross-validation to get coefficient set:
    for j = 1:k
        % split up CV data sets:
        these_CV_indices = CV_indices(j,:);
        these_CV_indices(isnan(these_CV_indices)) = [];
        CV_valid_y = y(these_CV_indices); % CV validation pigments
        CV_train_y = y;
        CV_train_y(these_CV_indices) = []; % CV training pigments
        CV_valid_x = x(these_CV_indices,:); % CV validation spectra
        CV_train_x = x;
        CV_train_x(these_CV_indices,:) = []; % CV training spectra

        % Standardize spectra for pc's:
        CV_train_x_z = zscore(CV_train_x);

        % take principal components of this CV training set:
        [CV_EOFs_train,CV_AFs_train,~,~,~] = pca(CV_train_x_z,'Centered',false,'NumComponents',max_components,'Rows','complete');
        
        % find the best possible model in terms of # components to use for
        % this CV:
        for l = 1:size(CV_AFs_train,2)
            the_lin_mdl = fitlm(CV_AFs_train(:,1:l), CV_train_y); % for an MLR model w/ derivative components
    
            these_betas = the_lin_mdl.Coefficients; % pull out linear coefficients of AF's
            these_betas = table2array(these_betas(:,1)); % convert data table to a matlab array
            this_alpha = these_betas(1); % pull out intercept
            these_betas(1) = []; % remove intercept so you only have multipliers in these_betas
    
            % zscored spectral derivative space:
            these_betas = CV_EOFs_train(:,1:l) * these_betas;

            % un-zscore:
            this_alpha = this_alpha - sum(these_betas .* (mean(CV_train_x)' ./ std(CV_train_x)'));
            these_betas = these_betas ./ std(CV_train_x)';    

            % Model Validation:
            CV_modeled_pigs = (CV_valid_x * these_betas) + this_alpha; 
    
            % constrain modeled pigments based on the pft_index you're
            % modeling
            CV_modeled_pigs(CV_modeled_pigs < 0) = 0;
    
            percent_errors(1:length(CV_valid_y),l) = ((CV_valid_y - CV_modeled_pigs)./CV_valid_y).*100;
            mean_percent_error(l) = mean(abs(percent_errors(:,l)));
            median_percent_error(l) = median(abs(percent_errors(:,l)));
    
            % fit linear model to look at modeled vs. observed:
            lin_mdl_4_validation = fitlm(CV_valid_y,CV_modeled_pigs);
            R2s(l) = lin_mdl_4_validation.Rsquared.Ordinary;
            RMSEs(l) = lin_mdl_4_validation.RMSE;
        end
    
        if strcmp(mdl_pick_metric, 'R2') == 1
            n_modes_to_use(j) = find(R2s == max(R2s)); % for selecting based on R^2 of predicted vs. observed
        elseif strcmp(mdl_pick_metric, 'RMSE') == 1
            n_modes_to_use(j) = find(RMSEs == min(RMSEs)); % for selecting based on minimizing RMSE
        elseif strcmp(mdl_pick_metric, 'avg') == 1
            n_modes_to_use(j) = find(mean_percent_error == min(mean_percent_error)); % for selecting based on minimizing avg percent error
        elseif strcmp(mdl_pick_metric, 'med') == 1
            n_modes_to_use(j) = find(median_percent_error == min(median_percent_error)); % for selecting based on minimizing median percent error
        end
    
        % apply your optimized model and record the g.o.f. statistics for this k-th CV:
        the_lin_mdl = fitlm(CV_AFs_train(:,1:n_modes_to_use(j)), CV_train_y); % for an MLR model w/ derivative components
        these_betas = the_lin_mdl.Coefficients;
        these_betas = table2array(these_betas(:,1));
        this_alpha = these_betas(1); %intercept
        these_betas(1) = []; % model coefficients for each spectral AF

        % now turn model coefficients for AF's into coefficients for the combined
        % derivative spectra:
        these_betas = CV_EOFs_train(:,1:n_modes_to_use(j)) * these_betas;

        % un-zscore:
        alpha(j) = this_alpha - sum(these_betas .* (mean(CV_train_x)' ./ std(CV_train_x)'));
        betas(:,j) = these_betas ./ std(CV_train_x)';    
    end

    betas = mean(betas, 2);
    alpha = mean(alpha);
        
end

end
