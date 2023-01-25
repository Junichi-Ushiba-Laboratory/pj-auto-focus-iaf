function data_EEG = BayesIAFest(data_EEG,samp_target)
if nargin < 2
    samp_target = 60;
    
end
fprintf('start BayesIAFest %d %d %d\n', data_EEG.expcond);
list_var_target  = {'onlineIAF_SNR_wRem','onlineIAF_SNR_woRem'};
list_var_estimate= {'IAF_B_est','IAF_B_est_raw'};
num_var          = numel(list_var_target);
data_EEG.res_est = [];
x                = data_EEG.FOI_IAF;

for i_var = 1 : num_var
    if data_EEG.flag_band == 1
        gb        = Gauss1D_Bayes;
    elseif data_EEG.flag_band == 2
        gb        = Gauss1D_Beta;
    end
    data_samp = data_EEG.(list_var_target{i_var});
    data_samp = data_samp(:,1);
    sec_rest  = numel(data_samp);
    flag_est  = 0;
    IAFlist   = zeros(5,1);
    for i_samp = 1 : sec_rest
        FOI_samp = data_samp(i_samp);
        gb       = gb.learn(FOI_samp);
        gb        = gb.calc_list_sigmas;
        pdf1      = gb.getPDF(i_samp,x);
        [IAF,prob]= gb.decideIAF(pdf1,x);
        IAFlist   = [IAFlist(2:end);IAF];
        if numel(unique(IAFlist)) == 1 && flag_est == 0 && i_samp >= samp_target
            data_EEG.(list_var_estimate{i_var}) = [IAF,prob,i_samp]';
            flag_est  = 1;
        end
    end
    gb               = gb.calc_list_sigmas;
    data_EEG.res_est = [data_EEG.res_est,gb];
    data_EEG         = data_EEG.decideIAF_B(gb,i_var);
end
fprintf('finish BayesIAFest %d %d %d\n', data_EEG.expcond);
end