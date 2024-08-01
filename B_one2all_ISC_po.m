clc;
clear all;

dirs.tlrc = '../../data/transformed_mat/positive12dofsm4_3mm_t/tlrc_t/';
dirs.positiveresult = '../../T_RESULTS/betweengroup/positiveresult/ISC';
savepath = fullfile(dirs.positiveresult);
if ~exist(savepath), mkdir(savepath); end

% addpath(genpath('../9_NIFTI_tools'));
addpath(genpath('../help_scripts'));

load(fullfile(dirs.tlrc, 'ALL_tlrc_po.mat'));

all_tlrc = vertcat(all_tlrc1,all_tlrc2);

All_To_other_ISC = [];
sub_id = [];

% Everyone to allother ISC
for i = 1:60
    
    fprintf('subject %i \n', i);
    subj1 = sprintf('t%03d', i);
    subj1_2other_ISC = [];
    
        for j = 1:60
        subj2 = sprintf('t%03d', j);
        if strcmp(subj1, subj2)
            continue;
        end
        target_subject = cellstr(strcat(subj1, '_', subj2));
        kept_id = (strcmp(tlrc_id(:),target_subject));
        subj1_subj2_r = all_tlrc((kept_id),1:end);
        subj1_2other_ISC = vertcat(subj1_2other_ISC, subj1_subj2_r);
        end
    subj1_2other_ISC = nanmean(subj1_2other_ISC);
    
    sub_id{i} = subj1;
    All_To_other_ISC(i,:) = subj1_2other_ISC;
end

All_To_other_ISC = All_To_other_ISC';
%calculate the average one2allother ISC
Ave_one2allother_ISC = nanmean(All_To_other_ISC,2);

save(fullfile(dirs.positiveresult, 'One2allother_ISC_po.mat'),'sub_id','All_To_other_ISC','Ave_one2allother_ISC');

%% Shuffle stats 
load(fullfile(dirs.positiveresult,'One2allother_ISC_po.mat'));
true_mean = Ave_one2allother_ISC;  % 

r_count = zeros(length(Ave_one2allother_ISC), 1);

for iteration = 1:10000
    
    fprintf('Iteration %i \n', iteration);
    
    % Random sign flipping 
    sign_flip = (rand(size(All_To_other_ISC)) > 0.5) * 2 - 1;
    fake_concat_corr = sign_flip .* All_To_other_ISC;
    fake_mean = nanmean(fake_concat_corr,2);  
    
    this_count = fake_mean > true_mean;
    r_count = r_count + this_count; 
    
end

p_permutation = (r_count + 1)/ iteration;
p_permutation(p_permutation > 1) = 1.00;

save(fullfile(dirs.positiveresult,'One2allother_ISC_po_flip10000.mat'),'true_mean','r_count','iteration','p_permutation');  

%% Save avg_one2allother map 
load(fullfile(dirs.positiveresult,'One2allother_ISC_po_flip10000.mat'));

% Create data for avg_one2all map 
% nii = load_nii('../../data/masks/3mm_standard.nii');
nii = load_nii('s001_s002_tlrc.nii');
nii.hdr.dime.datatype = 64;
nii.hdr.dime.glmax = -1;

data = true_mean;
nii.img = data;
nii.img(isnan(nii.img)) = 0;

save_nii(nii, fullfile(dirs.positiveresult,'Avg_one2allother_Nocorrection.nii'));

%% Calculate p_value and run FDR correction 

[h, crit_p] = fdr_bky(p_permutation, 0.002, 'yes');

% Convert back to full space 
p_mask = h;
nii.img(~p_mask) = 0;   

save_nii(nii, fullfile(dirs.positiveresult,'Avg_one2allother_fdr0002.nii'));

[h, crit_p] = fdr_bky(p_permutation, 0.001, 'yes');

% Convert back to full space 
p_mask = h;
nii.img(~p_mask) = 0;   

save_nii(nii, fullfile(dirs.positiveresult,'Avg_one2allother_fdr0001.nii'));