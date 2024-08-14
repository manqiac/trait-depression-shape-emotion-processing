clc;
clear all;

dirs.tlrc = '../../data/transformed_mat/positive12dofsm4_3mm_t/tlrc_t/';
dirs.positiveresult = '../../T_RESULTS/betweengroup/positiveresult/ISC';
savepath = fullfile(dirs.positiveresult);
if ~exist(savepath), mkdir(savepath); end

% addpath(genpath('../9_NIFTI_tools'));
addpath(genpath('../help_scripts'));

load(fullfile(dirs.tlrc, 'ALL_tlrc_po.mat'));

WithinISC = [];
BetweenISC = [];
WithinBetweenISC = [];
sub = [];

% LT_withinISC & LT_betweenISC
for i = 1:30
    fprintf('subject %i \n', i);
    
    subj1 = sprintf('t%03d', i);
    subj1_withinISC = [];
    subj1_betweenISC = [];
        for j = 1:30
            subj2 = sprintf('t%03d', j);
            % if subj1 == subj2
            if strcmp(subj1, subj2)
                continue;
            end
            target_subject = strcat(subj1, '_', subj2);        
            kept_id = (strcmp(tlrc_id(:),target_subject));        
            subj1_subj2_r = all_tlrc((kept_id),1:end);        
            subj1_withinISC = vertcat(subj1_withinISC, subj1_subj2_r);
        end
        
        subj1_meanwithinISC = nanmean(subj1_withinISC);
        
        for j = 31:60
            subj2 = sprintf('t%03d', j);
            target_subject = strcat(subj1, '_', subj2);
            kept_id = (strcmp(tlrc_id(:),target_subject));
            subj1_subj2_r = all_tlrc((kept_id),1:end);
            subj1_betweenISC = vertcat(subj1_betweenISC, subj1_subj2_r);
        end
        
        subj1_meanbetweenISC = nanmean(subj1_betweenISC);
        subj1_meanwithinbetweenISC = subj1_meanwithinISC - subj1_meanbetweenISC;
        
        sub{i} = subj1;
        WithinISC(i,:) = subj1_meanwithinISC;
        BetweenISC(i,:) = subj1_meanbetweenISC;
        WithinBetweenISC(i,:) = subj1_meanwithinbetweenISC;
end

% HT_withinISC & HT_betweenISC
for i = 31:60
    subj1 = sprintf('t%03d', i);
    fprintf('subject %i \n', i);
    
    subj1_withinISC = [];
    subj1_betweenISC = [];
        for j = 1:30
        subj2 = sprintf('t%03d', j);
        target_subject = strcat(subj1, '_', subj2);
        kept_id = (strcmp(tlrc_id(:),target_subject));%
        subj1_subj2_r = all_tlrc((kept_id),1:end);
        subj1_betweenISC = vertcat(subj1_betweenISC, subj1_subj2_r);
        end
        
        subj1_meanbetweenISC = nanmean(subj1_betweenISC);
        
        for j = 31:60
            subj2 = sprintf('t%03d', j);
            if strcmp(subj1, subj2)           
                continue;
            end
            target_subject = strcat(subj1, '_', subj2);
            kept_id = (strcmp(tlrc_id(:),target_subject));
            subj1_subj2_r = all_tlrc((kept_id),1:end);
            subj1_withinISC = vertcat(subj1_withinISC, subj1_subj2_r);
        end
        
        subj1_meanwithinISC = nanmean(subj1_withinISC);
        subj1_meanwithinbetweenISC = subj1_meanwithinISC - subj1_meanbetweenISC;
        
        sub{i} = subj1;
        WithinISC(i,:) = subj1_meanwithinISC;
        BetweenISC(i,:) = subj1_meanbetweenISC;
        WithinBetweenISC(i,:) = subj1_meanwithinbetweenISC;
end

save(fullfile(dirs.positiveresult, 'WithinBetweenISC_po.mat'),'sub','WithinISC','BetweenISC','WithinBetweenISC');

%% Shuffle stats 
load(fullfile(dirs.positiveresult,'WithinBetweenISC_po.mat'));

% brain mask
mask_struc = load_nii('s001_s002_tlrc.nii');
mask_struc.img(isnan(mask_struc.img)) = 0;
mask_4d = mask_struc.img;
mask_dimensions = size(mask_struc.img);
mask = logical(reshape(mask_4d, [mask_dimensions(1) * mask_dimensions(2) * mask_dimensions(3), 1]));

%% within-between ISC
WithinBetweenISC = WithinBetweenISC';
within_between_mask = zeros(length(WithinBetweenISC(:,1)), 60);
within_between_mask(~mask,:) = 0;
within_between_mask(mask,:) = WithinBetweenISC(mask,:);
    
% Shuffle stats 
true_mean = nanmean(within_between_mask,2);  % 
r_count = zeros(length(true_mean), 1);
    
for iteration = 1:10000
    fprintf('Iteration %i \n', iteration);
        
    % Random sign flipping 
    sign_flip = (rand(size(within_between_mask)) > 0.5) * 2 - 1;
    fake_within_between_mask = sign_flip .* within_between_mask;
    fake_mean = nanmean(fake_within_between_mask,2);  ?
        
    this_count = fake_mean > true_mean;
    r_count = r_count + this_count; 
end
    
    permutation_p = (r_count + 1) / iteration;
    permutation_p(permutation_p > 1) = 1.0000;
   
    save(fullfile(dirs.positiveresult,sprintf('WithinBetweenISC_flip10000_po.mat')),'true_mean','permutation_p','r_count','iteration');
    

%% Save map 
load(fullfile(dirs.positiveresult,'WithinBetweenISC_flip10000_po.mat'));

% Create data for avg_one2all map 
nii = load_nii('s001_s002_tlrc.nii');
nii.hdr.dime.datatype = 64;
nii.hdr.dime.glmax = -1;

data = true_mean;
nii.img = data;
nii.img(isnan(nii.img)) = 0;

save_nii(nii,fullfile(dirs.positiveresult,'WithinBetweenISC_R_Nocorrection.nii'));

%% Calculate p_value and run FDR correction 

[h, crit_p] = fdr_bky(permutation_p, 0.002, 'yes');

% Convert back to full space 
p_mask = h;
nii.img(~p_mask) = 0;   

save_nii(nii, fullfile(dirs.positiveresult,'WithinBetweenISC_R_fdr0002.nii'));

% Save z-map 
z = @(p) -sqrt(2) * erfcinv(p*2);
permutation_p(permutation_p == 0) = 0.0001;
z_value = z(1-permutation_p);
z_data = z_value;
nii.img = z_data;
nii.img(isnan(nii.img)) = 0;
nii.img(~p_mask) = 0;

save_nii(nii,fullfile(dirs.positiveresult,'WithinBetweenISC_Z_fdr0002.nii'));

%% Calculate p_value and run FDR correction 

[h, crit_p] = fdr_bky(permutation_p, 0.001, 'yes');

% Convert back to full space 
p_mask = h;
nii.img(~p_mask) = 0;   

save_nii(nii, fullfile(dirs.positiveresult,'WithinBetweenISC_R_fdr0001.nii'));

% Save z-map
z = @(p) -sqrt(2) * erfcinv(p*2);
permutation_p(permutation_p == 0) = 0.0001;
z_value = z(1-permutation_p);
z_data = z_value;
    
nii.img = z_data;
nii.img(isnan(nii.img)) = 0;
nii.img(~p_mask) = 0;

save_nii(nii,fullfile(dirs.positiveresult,'WithinBetweenISC_Z_fdr0001.nii'));

