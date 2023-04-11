function [smooth_spect] = smooth_spectra(spectra, filt_type, filt_width)
% applies smoothing filter to spectra
% assumes spectra are sampled uniformly
% input spectra is an observation x spectra matrix. wl is a vector of
% length size(spectra,2)
% filt_type is one of 'moving', 'sgolay' (Savitsky-Golay), 'rlowess',
% 'rloess', 'hamming
% filt_width is a scalar or vector or filter widths (in units of wl)
% output smooth_spec can assume 2 types: if filt_width is scalar, it is a
% matrix matching spectra; if filt_width is vector, it is a cell array of
% matrices with smoothed spectra in each element using the corresponding
% value of filt_width

all_smooth = cell(length(filt_width), 1);
for i = 1:length(filt_width)

    fw = filt_width(i);
    if strcmp(filt_type, 'hamming')
        % hamming is special ... not supported by built-in smooth
        % function
        % design hamming window filter:
        d = designfilt('lowpassfir','FilterOrder',fw-1,...
            'CutoffFrequency',(1/fw),'DesignMethod','window',...
            'Window','hamming','SampleRate',1); 

        b = d.Coefficients; % filter coefficients for convolution

        smooth_spect = nan(size(spectra)); % initialize smoothed matrix
        for j = 1:size(spectra, 1)
            smooth_spect(j, :) = conv(spectra(j, :), b, 'same');
        end
    else
        % use built-in smooth function:
        smooth_spect = nan(size(spectra)); % initialize smoothed matrix
        for j = 1:size(spectra, 1)
            smooth_spect(j, :) = smooth(spectra(j, :), fw, filt_type);
        end
    end

    % fill in unfiltered edges with NaNs
    smooth_spect(: , 1:((fw/2)-0.5)) = NaN; 
    smooth_spect(: , end-((fw/2)-1.5):end) = NaN; 
    all_smooth{i} = smooth_spect;
end

if length(all_smooth) == 1
    % de-cell array and return a matrix:
    smooth_spect = all_smooth{1};
else
    % just change the name:
    smooth_spect = all_smooth;
end


