clc;
clear all;

dirs.mat = '../../T_RESULTS/betweengroup/positiveresult/ISC';
dirs.mask = '../../T_RESULTS/betweengroup/positiveresult/ISC_maskforISFC/';
roi_wbISC = [];
rois = {'Temporal_Inf_Lmask' 'Amygdala_Lmask' 'Insula_Rmask' 'Frontal_Inf_Orb_Lmask' 'Frontal_Inf_Orb_Rmask' 'Frontal_Sup_Lmask' 'Cingulum_Mid_Rmask' 'Precuneus_Rmask' 'SupraMarginal_Lmask' 'Paracentral_Lobule_Lmask' 'Frontal_Mid_Lmask' 'Parietal_Sup_Rmask' 'Postcentral_Lmask' 'Supp_Motor_Area_Rmask' 'Precentral_Rmask' 'Precentral_Lmask' 'Frontal_Sup_Medial_Rmask' 'Precuneus_Rmask'};

for i =1:length(rois)
    roi=rois{i};
    
    load(fullfile(dirs.mat,'WithinBetweenISC_po.mat'));
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
csvwrite(fullfile(dirs.mask, 'WithinBetweenISC_ROI.csv'),roi_wbISC);