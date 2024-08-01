% Before this script,the masks/standard and transformed_data should be
% prepared. And you should install Add-Ons Tools for NIFTI and ANALYZE image.
%   Get data into .mat files for faster analysis

clc;
clear all;

dirs.fMRI = '../../data/processed_data/Nega_voxel_processed';   % 
dirs.mat = '../../data/transformed_mat/negative12dofsm4_3mm_t';   %

if ~exist(dirs.mat), mkdir(dirs.mat); end;

addpath(genpath('../help_scripts'));

%% transfer movie2.nii to movie2.mat file

for i = 1:60
    subj1 = sprintf('t%03d', i);
    fprintf('Running Subject %i \n', i);
    
    
    nii = load_nii(fullfile(dirs.fMRI,sprintf('%s.nii', subj1)));        
    data = nii.img;         
    datasize = size(data);        
        
    data = single(reshape(data, [(size(data, 1) * size(data, 2) * size(data, 3)), size(data, 4)]));   % ��������        

    A = data';
    movie_data = zscore(A);
    movie_data = movie_data';
      
    save(fullfile(dirs.mat, sprintf('%s.mat', subj1)), 'movie_data', 'data', 'datasize');
end

%% Extract ROI timecourse for each participant
dirs.roi = '../../T_RESULTS/betweengroup/negativeresult/ISC_maskforISFC/'; 
dirs.roi_tc = '../../data/transformed_mat/negative12dofsm4_3mm_t/roi_tc/';
roi_tcpath = fullfile(dirs.roi_tc);
if ~exist(roi_tcpath), mkdir(roi_tcpath); end;

rois = {'Temporal_Inf_Lmask' 'Hippocampus_Rmask' 'Lingual_Rmask' 'Putamen_Rmask' 'SupraMarginal_Rmask' 'Precuneus_Rmas' 'Precuneus_R2mask' 'Frontal_Sup_Lmask' 'Precentral_Lmask'};

for r = 1:length(rois)
    roi = rois{r};  
    roi_tc = [];
    % Load mask
    mask_struc = load_nii(fullfile(dirs.roi,sprintf('%s.nii',roi)));
    mask_4d = mask_struc.img;
    mask_dimensions = size(mask_struc.img);
    mask = logical(reshape(mask_4d,[mask_dimensions(1)*mask_dimensions(2)*mask_dimensions(3),1]));

    for i = 1:60
         subj1 = sprintf('t%03d', i);
         fprintf('Running Subject %i \n', i);
        % Load subject's data
    
        load(fullfile(dirs.mat, sprintf('%s.mat', subj1)));
        % Mask data
        B = movie_data(mask,:);
        roi_tc(:,i) = mean(B);

    end
    save(fullfile(dirs.roi_tc,sprintf('%s.mat', roi)), 'roi_tc');
end

%% ISFC analysis
dirs.result = '../../T_RESULTS/betweengroup/negativeresult/ISFC/';
if ~exist(fullfile(dirs.result)), mkdir(fullfile(dirs.result)); end;

%% Fisher r-to-z
All_within_ISFC = [];
All_between_ISFC = [];

ROI_pair = cell(1, nchoosek(length(rois), 2));
pairIndex = 1;

for r1 = 1:length(rois)
            roi1 = rois{r1};
            load(fullfile(dirs.roi_tc,sprintf('%s.mat', roi1)));
            roi1_tc = roi_tc;
            between_ISFC = [];
            within_ISFC = [];
            for r2 = r1+1:length(rois);
                roi2 = rois{r2};
                load(fullfile(dirs.roi_tc,sprintf('%s.mat', roi2)));
                roi2_tc = roi_tc;
                for s1 = 1:30               
                    subj1 = sprintf('t%03d', s1);                
                    fprintf('Running subject %s\n',subj1);
                    subj1_between_ISFC = [];                
                    subj1_within_ISFC = [];                              
                    
                    for s2 = 1:30                    
                        subj2 = sprintf('t%03d', s2);                    
                        if strcmp(subj1, subj2)                     
                            continue;
                        end
                        subj1_within_ISFC(s1,s2) = corr(roi1_tc(:,s1), roi2_tc(:,s2));     
                    end
                    subj1_within_ISFC(find(subj1_within_ISFC==0))=[];
                    within_ISFC(s1,r2) = mean(subj1_within_ISFC,2);
                    for s2 = 31:60
                        subj2 = sprintf('t%03d', s2);                    
                        subj1_between_ISFC(s1,s2) = corr(roi1_tc(:,s1), roi2_tc(:,s2));                                        
                    end
                    subj1_between_ISFC(find(subj1_between_ISFC==0))=[];
                    between_ISFC(s1,r2) = mean(subj1_between_ISFC,2);   
                end

                for s1 = 31:60
                    subj1 = sprintf('t%03d', s1);                
                    fprintf('Running subject %s \n', subj1);
                    subj1_between_ISFC = [];                
                    subj1_within_ISFC = [];                              
                    for s2 = 1:30                    
                        subj2 = sprintf('t%03d', s2);                                                                      
                        subj1_between_ISFC(s1,s2) = corr(roi1_tc(:,s1), roi2_tc(:,s2));     
                    end
                    subj1_between_ISFC(find(subj1_between_ISFC==0))=[];
                    between_ISFC(s1,r2) = mean(subj1_between_ISFC,2);
                    for s2 = 31:60
                        subj2 = sprintf('t%03d', s2);                                                
                        if strcmp(subj1, subj2)                     
                            continue;
                        end      
                        subj1_within_ISFC(s1,s2) = corr(roi1_tc(:,s1), roi2_tc(:,s2));                                            
                    end
                    subj1_within_ISFC(find(subj1_within_ISFC==0))=[];
                    within_ISFC(s1,r2) = mean(subj1_within_ISFC,2);
                end
                roicombined = strcat(roi1, '_', roi2);        
                ROI_pair{pairIndex} = roicombined;        
                pairIndex = pairIndex + 1;                            
            end
            All_within_ISFC = horzcat(All_within_ISFC, within_ISFC);
            All_between_ISFC = horzcat(All_between_ISFC, between_ISFC);          
end

All_within_ISFC(:,all(All_within_ISFC == 0)) =[];
All_between_ISFC(:,all(All_between_ISFC == 0)) =[];
All_within_ISFC = atanh(All_within_ISFC);             
All_between_ISFC = atanh(All_between_ISFC);
All_withinbetween_ISFC  = All_within_ISFC - All_between_ISFC;

save(fullfile(dirs.result,sprintf('WithinBetweenISFC_nega.mat')),'ROI_pair','All_within_ISFC','All_between_ISFC','All_withinbetween_ISFC');
csvwrite(fullfile(dirs.result, 'All_within_ISFC_ROIwise_nega.csv'),All_within_ISFC);
csvwrite(fullfile(dirs.result, 'All_between_ISFC_ROIwise_nega.csv'),All_between_ISFC);
csvwrite(fullfile(dirs.result, 'All_withinbetween_ISFC_ROIwise_nega.csv'),All_withinbetween_ISFC);