clc;
clear all;

dirs.mat = '../../T_RESULTS/betweengroup/negativeresult/ISC';
dirs.mask = '../../T_RESULTS/betweengroup/negativeresult/ISC_maskforISFC/';
roi_wbISC = [];
rois = {'Temporal_Inf_Lmask' 'Hippocampus_Rmask' 'Lingual_Rmask' 'Putamen_Rmask' 'SupraMarginal_Rmask' 'Precuneus_Rmas' 'Precuneus_R2mask' 'Frontal_Sup_Lmask' 'Precentral_Lmask'};

for i =1:length(rois)
    roi=rois{i};
    
    load(fullfile(dirs.mat,'WithinBetweenISC_nega.mat'));
    % Load mask
    mask_struc = load_nii(fullfile(dirs.mask,sprintf('%s.nii',roi)));
    mask_4d = mask_struc.img;
    mask_dimensions = size(mask_struc.img);
    mask = logical(reshape(mask_4d,[mask_dimensions(1)*mask_dimensions(2)*mask_dimensions(3),1]));
    
    % Mask data
    WithinBetweenISC = WithinBetweenISC';
    A = WithinBetweenISC(mask,:);
    roi_wbISC(i,:) = mean(A);
    
end

roi_wbISC = roi_wbISC';

save(fullfile(dirs.mask,'WithinBetweenISC_ROI.mat'),'rois','roi_wbISC');
csvwrite(fullfile(dirs.mask, 'WithinBetweenISC_ROI.csv'),roi_wbISC)
