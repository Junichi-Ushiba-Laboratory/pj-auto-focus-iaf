function [out,fit] = fcn_removeF(in,fint)
    if nargin < 2
        fint = [8 13];
    else
        fint = [fint(1) fint(end)];
    end
    frequency_interval = fint;
    fmax               = 50; % 65
    power              = in;
    frequency          = [1 : numel(power)]';
    bins_selected      = find(frequency>=3,1):find(frequency>=45,1);
    
    frequency     = frequency(bins_selected)';
    power         = power(bins_selected,:);
    
    log_freq      = log10(frequency);
    log_power     = 10*log10(power);
    
    % define frequence bands
    idx1 = find(frequency>=0.5 & frequency <= 7);
    idx2 = find(frequency>=35 & frequency <= fmax);
    idx3 = find(frequency>=frequency_interval(1) & frequency <= frequency_interval(2));
    idx4 = find(frequency>=.5 & frequency<fmax);
    
    idx_freq        = [idx1 idx2];
    log_freq_roi    = [ones(1,numel(idx_freq)); log_freq(idx_freq)];
    log_power_roi   = [log_power(idx1)' log_power(idx2)'];
    
    % fit linear slope in log-log space
    
    fit         = log_freq_roi'\log_power_roi';
    slope       = fit(2);
    intercept   = fit(1);
    fit_1f      = slope*log10(frequency') + intercept;
    out         = in;
    out(bins_selected) = log_power - fit_1f;
end