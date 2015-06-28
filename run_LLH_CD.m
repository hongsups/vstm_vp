function run_LLH_CD(table_ID,capacity_idx,expid,subjid)

% compute LLH and ML estimate                

file_path = 'CD/summary/';
sbj_name = [num2str(expid) 'S' num2str(subjid)];
[~, model_name] = create_modelname([],[],[],[],capacity_idx);
if table_ID == 1
    create_LLH_IL_CD(capacity_idx,expid,subjid);
    close all
elseif table_ID == 2
    fname = [file_path 'LLH_SA_' model_name '_' sbj_name '.mat'];
%     if ~exist(fname,'file')
        [L, ML_estimate, BIC] = create_LLH_SA_CD(capacity_idx,expid,subjid);
        save(fname,'L','ML_estimate','BIC');  
%     else
%         fprintf('File already exists - skipping\n');
%     end
elseif table_ID == 3
    fname = [file_path 'LLH_EP_' model_name '_' sbj_name '.mat'];
%     if ~exist(fname,'file')
        [L, ML_estimate, BIC] = create_LLH_EP_CD(capacity_idx,expid,subjid);
        save(fname,'L','ML_estimate','BIC');  
%     else
%         fprintf('File already exists - skipping\n');
%     end
else
    fname = [file_path 'LLH_VP_' model_name '_' sbj_name '.mat'];
%     if ~exist(fname,'file')
        [L, ML_estimate, BIC] = create_LLH_VP_CD(capacity_idx,expid,subjid);
        save(fname,'L','ML_estimate','BIC');  
%     else
%         fprintf('File already exists - skipping\n');
%     end
end

