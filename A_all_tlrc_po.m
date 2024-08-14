clc;
clear all;

dirs.fMRI = '../../data/processed_data/Posi_voxel_processed/tlrc_t/';
dirs.tlrc = '../../data/transformed_mat/positive12dofsm4_3mm_t/tlrc_t/';
tlrc_outputpath = fullfile(dirs.tlrc);
if ~exist(tlrc_outputpath), mkdir(tlrc_outputpath); end


all_tlrc = [];
tlrc_id = [];

for i = 1:60
    subj1 = sprintf('t%03d', i);
    fprintf('subject %i \n', i);
    for j = i+1:60
        subj2 = sprintf('t%03d', j);
        nii = load_nii(fullfile(dirs.fMRI,sprintf('%s_%s_tlrc.nii', subj1, subj2)));  
        if strcmp(subj1, subj2)
            continue;
        end
        data = nii.img;  
        datasize = size(data);
        data = single(reshape(data, [(size(data, 1) * size(data, 2) * size(data, 3)), size(data, 4)]));   % 重塑数据
%       keptvox = mdata > mcutoff;   
%       movie_data = data(keptvox, :);
        subj1_subj2 = cellstr(strcat(subj1, '_', subj2));
        subj2_subj1 = cellstr(strcat(subj2, '_', subj1));
        a = [subj1_subj2; subj2_subj1];
        b = vertcat(data',data');
        tlrc_id = vertcat(tlrc_id, a);
        all_tlrc = vertcat(all_tlrc, b);
    end
end

all_tlrc = all_tlrc(1:2:end,:);

save(fullfile(dirs.tlrc, 'ALL_tlrc_po.mat'),'tlrc_id','all_tlrc');
