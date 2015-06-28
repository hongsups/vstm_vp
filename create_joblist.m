clear all

% % jobllist for lookup table
% expid = input('EXP ID (orientation(1), color(2))..............:');
% tableid = input('MODEL ID (SA(1), EP(2), VP(3))..............:');
% 
% noise_idx = 1;
% capacity_idx = 1;
% param_info;
% 
% if tableid == 1
%     fid = fopen('joblist_SA.txt','w');    
%     for J1_SA_idx = 1:nsteps
%         for noise_idx = 1
%             fprintf(fid,'null,%d,%d,%d\n',J1_SA_idx,noise_idx,expid);
%         end
%     end
% elseif tableid == 2
%     fid = fopen('joblist_EP.txt','w');
%     for J1_EP_idx = 1:nsteps
%         for noise_idx = 1
%             fprintf(fid,'null,%d,%d,%d\n',J1_EP_idx,noise_idx,expid);
%         end
%     end
% elseif tableid == 3
%     fid = fopen('joblist_VP.txt','w');
%     for J1_VP_idx = 1:nsteps
%         for vari_idx = 1
%             for rule_idx = 1
%                 fprintf(fid,'null,%d,%d,%d,2,1,1\n',J1_VP_idx,vari_idx,rule_idx);
%             end
%         end
%     end
% else
%     fid = fopen('joblist_VP.txt','w');
%     for J1_VP_idx = 1:nsteps
%         for vari_idx = 1
%             for rule_idx = 1
% %                 for kappa_idx = 1:2
% %                     for noise_idx = 1:2
%                         fprintf(fid,'null,%d,%d,%d,%d,%d,%d\n',J1_VP_idx,vari_idx,rule_idx,2,1,expid);
% %                     end
% %                 end
%             end
%         end
%     end
% end
% 
% fclose(fid);

%% create model name 
% 
% n1 = ['T' 'V'];
% n2 = ['O' 'S'];
% n3 = ['E' 'I'];
% n4 = ['N' 'X'];
% n5 = ['F' 'P' 'U' 'X'];
% 
% fid = fopen('VPmodellist.txt','w');
% for i = 1:2
%     for j = 1:2
%         for k = 1:2
%             for l = 1:2
%                 for m = 1:4
%                     fprintf(fid,'%c%c%c%c%c\n',n1(i),n2(j),n3(k),n4(l),n5(m));
%                 end
%             end
%         end
%     end
% end
% fid(close);
% 
% fid = fopen('EPmodellist.txt','w');
% for l = 1:2
%     for m = 1:4
%         fprintf(fid,'%c%c\n',n4(l),n5(m));
%     end
% end
% fid(close);
% 
% fid = fopen('SAmodellist.txt','w');
% for l = 1:2
%     for m = 1:3   
%         fprintf(fid,'%c%c\n',n4(l),n5(m));
%     end
% end
% fid(close);
% 
%% create job list for fake data test
% fid = fopen('joblist_all.txt','w');    
% for model_idx = 1:4
%     if model_idx == 1
%         for capacity_idx = 1:3
%             fprintf(fid,'null,%d,%d,%d,%d,%d,%d\n',model_idx,1,1,1,1,capacity_idx+1);
%         end
%     elseif model_idx == 2
%         for capacity_idx = 1:3
%             for noise_idx = 1:2
%                 fprintf(fid,'null,%d,%d,%d,%d,%d,%d\n',model_idx,1,1,1,noise_idx,capacity_idx+1);
%             end
%         end
%     elseif model_idx == 3
%         for capacity_idx = 1:4
%             for noise_idx = 1:2
%                 fprintf(fid,'null,%d,%d,%d,%d,%d,%d\n',model_idx,1,1,1,noise_idx,capacity_idx);
%             end
%         end
%     else
%         for vari_idx = 1:2
%             for rule_idx = 1:2
%                 for kappa_idx = 1:2
%                     for noise_idx = 1:2
%                         for capacity_idx = 1:4
%                             fprintf(fid,'null,%d,%d,%d,%d,%d,%d\n',model_idx,vari_idx,rule_idx,kappa_idx,noise_idx,capacity_idx);
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% fclose(fid);

% %% create fake data test
% subj_ids = {'QFA','QFB','QFC','QFD','QFE'};    
% for i = 1:length(subj_ids)
%     subjid = subj_ids{i};
%     for model_idx = 1:4
%         if model_idx == 1
%             for capacity_idx = 1:3
%                 [table_name model_name] = create_modelname([],[],[],[],capacity_idx+1);
%                 sbjname_IL{i,capacity_idx} = [num2str(model_idx) '_' model_name '_' subjid];
%             end
%         elseif model_idx == 2
%             for capacity_idx = 1:3
%                 for noise_idx = 1:2
%                     [table_name model_name] = create_modelname([],[],[],noise_idx,capacity_idx+1);
%                     sbjname_SA{i,noise_idx,capacity_idx} = [num2str(model_idx) '_' model_name '_' subjid];
%                 end
%             end
%         elseif model_idx == 3
%             for capacity_idx = 1:4
%                 for noise_idx = 1:2
%                     [table_name model_name] = create_modelname([],[],[],noise_idx,capacity_idx);
%                     sbjname_EP{i,noise_idx,capacity_idx} = [num2str(model_idx) '_' model_name '_' subjid];
%                 end
%             end
%         else
%             for vari_idx = 1:2
%                 for rule_idx = 1:2
%                     for kappa_idx = 1:2
%                         for noise_idx = 1:2
%                             for capacity_idx = 1:4
%                                 [table_name model_name] = create_modelname(vari_idx,rule_idx,kappa_idx,noise_idx,capacity_idx);
%                                 sbjname_VP{i,vari_idx,rule_idx,kappa_idx,noise_idx,capacity_idx} = [num2str(model_idx) '_' model_name '_' subjid];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% sbjname = [sbjname_IL(:); sbjname_SA(:); sbjname_EP(:); sbjname_VP(:)];
% fname = ('fake_subjects.mat');
% save(fname,'sbjname');

% %% create model name fake data test
% for model_idx = 1:4
%     if model_idx == 1
%         for capacity_idx = 1:3
%             [table_name model_name] = create_modelname([],[],[],[],capacity_idx+1);
%             sbjname_IL{capacity_idx} = [model_name];
%         end
%     elseif model_idx == 2
%         for capacity_idx = 1:3
%             for noise_idx = 1:2
%                 [table_name model_name] = create_modelname([],[],[],noise_idx,capacity_idx+1);
%                 sbjname_SA{noise_idx,capacity_idx} = [model_name];
%             end
%         end
%     elseif model_idx == 3
%         for capacity_idx = 1:4
%             for noise_idx = 1:2
%                 [table_name model_name] = create_modelname([],[],[],noise_idx,capacity_idx);
%                 sbjname_EP{noise_idx,capacity_idx} = [model_name];
%             end
%         end
%     else
%         for vari_idx = 1:2
%             for rule_idx = 1:2
%                 for kappa_idx = 1:2
%                     for noise_idx = 1:2
%                         for capacity_idx = 1:4
%                             [table_name model_name] = create_modelname(vari_idx,rule_idx,kappa_idx,noise_idx,capacity_idx);
%                             sbjname_VP{vari_idx,rule_idx,kappa_idx,noise_idx,capacity_idx} = [model_name];
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % models = [sbjname_IL(:); sbjname_SA(:); sbjname_EP(:); sbjname_VP(:)];
% sbjname_IL = sbjname_IL(:);
% sbjname_SA = sbjname_SA(:); 
% sbjname_EP = sbjname_EP(:);
% sbjname_VP = sbjname_VP(:);
% fname = ('modelnames_2.mat');
% save(fname,'sbjname_IL','sbjname_SA','sbjname_EP','sbjname_VP');

%% create job list for LLH VP (fake data test)
% fid = fopen('joblist_LLH_VP_fake.txt','w');    
%             for vari_idx = 1:2
%                 for rule_idx = 1:2
%                     for kappa_idx = 1:2
%                         for noise_idx = 1:2
%                             for capacity_idx = 1:4
%                                 for subj_idx = 1:405
%                                     fprintf(fid,'null,%d,%d,%d,%d,%d,%d,%d\n',vari_idx,rule_idx,kappa_idx,noise_idx,capacity_idx,3,subj_idx);
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
% fclose(fid);

% fid = fopen('joblist_LLH_EP_fake.txt','w');    
%     for noise_idx = 1:2
%         for capacity_idx = 1:4
%             for subj_idx = 1:405
%                 fprintf(fid,'null,%d,%d,%d,%d\n',noise_idx,capacity_idx,3,subj_idx);
%             end
%         end
%     end
% fclose(fid);

% fid = fopen('joblist_LLH_SA_fake.txt','w');    
%     for noise_idx = 1:2
%         for capacity_idx = 1:3
%             for subj_idx = 1:405
%                 fprintf(fid,'null,%d,%d,%d,%d\n',noise_idx,capacity_idx+1,3,subj_idx);
%             end
%         end
%     end
% fclose(fid);

% fid = fopen('joblist_LLH_IL_fake.txt','w');    
%     for capacity_idx = 1:3
%         for subj_idx = 1:405
%             fprintf(fid,'null,%d,%d,%d\n',capacity_idx+1,3,subj_idx);
%         end
%     end
% fclose(fid);


%%
% fid = fopen('joblist_cluster.txt','w');    
% for i = 1:50
%     fprintf(fid,'null,%d\n',i);
% end
% fclose(fid);

%% create model name fake data test

% for model_idx = 1:4
%     if model_idx == 1
%         for capacity_idx = 1:3
%             [table_name model_name] = create_modelname([],[],[],[],capacity_idx+1);
%             sbjname_IL{capacity_idx} = [num2str(model_idx) '_' model_name];
%         end
%     elseif model_idx == 2
%         for capacity_idx = 1:3
%             for noise_idx = 1
%                 [table_name model_name] = create_modelname([],[],[],noise_idx,capacity_idx+1);
%                 sbjname_SA{noise_idx,capacity_idx} = [num2str(model_idx) '_' model_name];
%             end
%         end
%     elseif model_idx == 3
%         for capacity_idx = 1:4
%             for noise_idx = 1
%                 [table_name model_name] = create_modelname([],[],[],noise_idx,capacity_idx);
%                 sbjname_EP{noise_idx,capacity_idx} = [num2str(model_idx) '_' model_name];
%             end
%         end
%     else
%         for vari_idx = 1
%             for rule_idx = 1
%                 for kappa_idx = 1
%                     for noise_idx = 1
%                         for capacity_idx = 1:4
%                             [table_name model_name] = create_modelname(vari_idx,rule_idx,kappa_idx,noise_idx,capacity_idx);
%                             sbjname_VP{vari_idx,rule_idx,kappa_idx,noise_idx,capacity_idx} = [num2str(model_idx) '_' model_name];
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % models = [sbjname_IL(:); sbjname_SA(:); sbjname_EP(:); sbjname_VP(:)];
% sbjname_IL = sbjname_IL(:);
% sbjname_SA = sbjname_SA(:); 
% sbjname_EP = sbjname_EP(:);
% sbjname_VP = sbjname_VP(:);
% sbjname = [sbjname_IL; sbjname_SA; sbjname_EP; sbjname_VP];
% fname = ('modelnames_14_fake.mat');
% save(fname,'sbjname_IL','sbjname_SA','sbjname_EP','sbjname_VP','sbjname');

%% joblist for check_model_fit.m
% %% create model name fake data test
% fid = fopen('joblist_check_model.txt','w');    
% for model_id = 1:4
%     if model_id < 3             % IL, SA
%         model_max_id = 3;   
%     elseif model_id == 3        % EP
%         model_max_id = 4;
%     else                        % VP
%         model_max_id = 16;      
%     end
%     for subid = 1:model_max_id      % # of models in a model groupo
%         for n_rep = 1:10            % # of fake subjects
%             fprintf(fid,'null,%d,%d,%d\n',model_id,subid,n_rep);
%         end
%     end
% end
% fclose(fid);
% 
%% joblist for run_fit
% 
% expid = input('exp (ori=1, color=2):');
% if expid == 1
%     fid = fopen('joblist_run_fit_ori.txt','w');
%     for i = 1:11
%         for model_idx = 4
%             fprintf(fid,'null,%d,%d,%d\n',model_idx,1,i);
%         end
%     end
%     fclose(fid);
% elseif expid == 2
%     fid = fopen('joblist_run_fit_col.txt','w');
%     for i = 1:9
%         for model_idx = 4
%             fprintf(fid,'null,%d,%d,%d\n',model_idx,2,i);
%         end
%     end
%     fclose(fid);
% end   

%% create model name for extended version of factorial
% for model_idx = 1:4
%     if model_idx == 1
%         for capacity_idx = 1:3
%             model_name = create_modelname_temp(model_idx,[],[],capacity_idx+1);
%             sbjname_IL{capacity_idx} = model_name;
%         end
%     elseif model_idx == 2
%         for capacity_idx = 1:3
%             model_name = create_modelname_temp(model_idx,[],[],capacity_idx+1);
%             sbjname_SA{capacity_idx} = model_name;
%         end
%     elseif model_idx == 3
%         for capacity_idx = 1:4
%             model_name = create_modelname_temp(model_idx,[],[],capacity_idx);
%             sbjname_EP{capacity_idx} = model_name;
%         end
%     else
%         for vari_idx = 1
%             for rule_idx = 1
%                 for capacity_idx = 1:4
%                     model_name = create_modelname_temp(model_idx,vari_idx,rule_idx,capacity_idx);
%                     sbjname_VP{vari_idx,rule_idx,capacity_idx} = model_name;
%                 end
%             end
%         end
%     end
% end
% 
% % models = [sbjname_IL(:); sbjname_SA(:); sbjname_EP(:); sbjname_VP(:)];
% sbjname_IL = sbjname_IL(:);
% sbjname_SA = sbjname_SA(:); 
% sbjname_EP = sbjname_EP(:);
% sbjname_VP = sbjname_VP(:);
% fname = ('modelnames_14_name.mat');
% save(fname,'sbjname_IL','sbjname_SA','sbjname_EP','sbjname_VP');

% %% Generate fake data for all models
% fid = fopen('joblist_ML_plots_SA.txt','w');    
% for model_id = 2
%     if model_id < 3             % IL, SA
%         capacity_idx = 2:4;   
%     else                        % EP, VP
%         capacity_idx = 1:4;
%     end
%     for ic = 1:length(capacity_idx)
%         for n_rep = 1:20            % # of fake subjects
%             fprintf(fid,'null,%d,%d,%d\n',model_id,capacity_idx(ic),n_rep);
%         end
%     end
% end
% fclose(fid);

% %% conjunction table
% fid = fopen('joblist_VP_conj_r1r2.txt','w');    
% nsteps = 15;
% for J_ori_bar_idx = 1:nsteps
%     for J_col_bar_idx = 1:nsteps
%         for tau1_bar_idx = 1:nsteps
%             fprintf(fid,'null,%d,%d,%d\n',J_ori_bar_idx,J_col_bar_idx,tau1_bar_idx);
%         end
%     end
% end
% fclose(fid);

%% joblist for conjunction fake data test

% fid = fopen('joblist_VPSR_fake.txt','w');    
% nsbj = 15;
% for i = 1:nsbj
%     fprintf(fid,'null,%d\n',i);
% end
% fclose(fid);

%% create SA table (FACTORIAL)
% fid = fopen('joblist_SA_new.txt','w');    
% nsteps = 20;
% for J1_SA_idx = 1:nsteps
%     fprintf(fid,'null,%d\n',J1_SA_idx);
% end
% fclose(fid);

%% create SA table (FACTORIAL)
% fid = fopen('joblist_SA_create_fakedata.txt','w');    
% nsbj = 15;
% capacity_idx = 1:3; % +1 for the real use (SA can't have capacity index of 1)
% for i = 1:length(capacity_idx)
%     for j = 1:nsbj
%         fprintf(fid,'null,%d,%d\n',i+1,j);
%     end
% end
% fclose(fid);

% create model name for extended version of factorial
% for model_idx = 1:4
%     if model_idx == 1
%         for capacity_idx = 1:3
%             model_name = create_modelname_temp(model_idx,capacity_idx+1);
%             sbjname_IL{capacity_idx} = model_name;
%         end
%     elseif model_idx == 2
%         for capacity_idx = 1:3
%             model_name = create_modelname_temp(model_idx,capacity_idx+1);
%             sbjname_SA{capacity_idx} = model_name;
%         end
%     elseif model_idx == 3
%         for capacity_idx = 1:4
%             model_name = create_modelname_temp(model_idx,capacity_idx);
%             sbjname_EP{capacity_idx} = model_name;
%         end
%     else
%         for capacity_idx = 1:4
%             model_name = create_modelname_temp(model_idx,capacity_idx);
%             sbjname_VP{capacity_idx} = model_name;
%         end
%     end
% end
% 
% % models = [sbjname_IL(:); sbjname_SA(:); sbjname_EP(:); sbjname_VP(:)];
% sbjname_IL = sbjname_IL(:);
% sbjname_SA = sbjname_SA(:); 
% sbjname_EP = sbjname_EP(:);
% sbjname_VP = sbjname_VP(:);
% sbjname = [sbjname_IL; sbjname_SA; sbjname_EP; sbjname_VP];
% 
% fname = ('modelnames_14_name_final.mat');
% save(fname,'sbjname_IL','sbjname_SA','sbjname_EP','sbjname_VP','sbjname');

%% run_param_recovery_test_conj
% % ga_conj_allmodels
% % create_fakedata_conj_allmodels
% % run_param_recovery_test_conj
% fid = fopen('run_param_recovery_test_conj','w');    
% ExpIdxVec = 1:3;
% nsbj = 10;
% for ExpIdx = 1:length(ExpIdxVec)
%     for SbjIdx = 1:nsbj
%         if ExpIdx == 1
%             fprintf(fid,'null,%d,%d,%d\n',ExpIdx,3,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d\n',ExpIdx,5,SbjIdx);
%         elseif ExpIdx == 2
%             fprintf(fid,'null,%d,%d,%d\n',ExpIdx,1,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d\n',ExpIdx,2,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d\n',ExpIdx,5,SbjIdx);
%         else
%             fprintf(fid,'null,%d,%d,%d\n',ExpIdx,4,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d\n',ExpIdx,5,SbjIdx);
%         end 
%     end
% end
% fclose(fid);

%% ga_conj_allmodels
% create_fakedata_conj_allmodels
% exp1: 8 subjects
% exp2: 8 subjects
% exp3: 5 subjects
% TestIdx,ExpIdx,ModelIdx,SbjIdx

% fid = fopen('joblist_ga_conj_allmodels.txt','w');    
% TestIdx = 0;
% ExpIdxVec = 1:3;
% nsbjVec = [8 8 5];
% for ExpIdx = 1:length(ExpIdxVec)
%     nsbj = nsbjVec(ExpIdx);
%     for SbjIdx = 1:nsbj
%         if ExpIdx == 1
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,3,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,5,SbjIdx);
%         elseif ExpIdx == 2
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,1,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,2,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,5,SbjIdx);
%         else
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,4,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,5,SbjIdx);
%         end 
%     end
% end
% fclose(fid);

% fid = fopen('joblist_ga_conj_allmodels.txt','w');    
% TestIdx = 0;
% ExpIdxVec = 1:3;
% nsbjVec = [8 8 5];
% for ExpIdx = 1:length(ExpIdxVec)
%     nsbj = nsbjVec(ExpIdx);
%     for SbjIdx = 1:nsbj
%         if ExpIdx == 1
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,3,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,5,SbjIdx);
%         elseif ExpIdx == 2
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,1,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,2,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,5,SbjIdx);
%         else
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,4,SbjIdx);
%             fprintf(fid,'null,%d,%d,%d,%d\n',TestIdx,ExpIdx,5,SbjIdx);
%         end 
%     end
% end
% fclose(fid);

% fid = fopen('joblist_ga_e1n8.txt','w');    
% TestIdx = 0;
% ExpIdxVec = 1;
% nsbjVec = 8;
% for ExpIdx = 1:length(ExpIdxVec)
%     nsbj = nsbjVec(ExpIdx);
%     for SbjIdx = 1:nsbj
%         PopSize = 64;
%         GenNum = 64;
%         SampleSize = 1280;
%         ntrials = 2400;
%         N = 8;
%         
%         fprintf(fid,'null,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',TestIdx,ExpIdx,3,SbjIdx,PopSize,GenNum,SampleSize,ntrials,N);
%         fprintf(fid,'null,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',TestIdx,ExpIdx,5,SbjIdx,PopSize,GenNum,SampleSize,ntrials,N);
%     end
% end
% fclose(fid);

% fid = fopen('joblist_ga_leak.txt','w');    
% TestIdx = 0;
% nsbj = 5;
% omega_idxVec = 0:2;
% power_idxVec = 0;
% popSize = 64;
% genNum = 64;
% sampleSize = 640;
% ntrials = 10000;
% 
% for SbjIdx = 1:nsbj
%     for i = 1:length(omega_idxVec)
%         for ii = 1:length(power_idxVec)
%             fprintf(fid,'null,%d,%d,%d,%d,%d,%d,%d,%d\n',TestIdx,SbjIdx,omega_idxVec(i),power_idxVec(ii),popSize,genNum,sampleSize,ntrials);
%         end
%     end
% end
% fclose(fid);

%% CD VP table

% fid = fopen('joblist_CD_VP.txt','w');    
% for i = 1:15
%     for j = 1:15
%         fprintf(fid,'null,%d,%d\n',i,j);
%     end
% end
% fclose(fid);

% %% de_delay
% fid = fopen('joblist_de.txt','w');    
% SbjIdxVec = 1:5;
% DelayIdxVec = [1 2 3 6];
% TestIdx = 0;
% for i = 1:length(SbjIdxVec)
%     SbjIdx = SbjIdxVec(i);
%     for j = 1:length(DelayIdxVec)
%         DelayIdx = DelayIdxVec(j);
%         fprintf(fid,'null,%d,%d,%d\n',SbjIdx,DelayIdx,TestIdx);
%     end
% end
% % fclose(fid);

%% de_delay
fid = fopen('joblist_de_fake.txt','w');    
SbjIdxVec = 1:50;
TestIdx = 0;
for i = 1:length(SbjIdxVec)
    fprintf(fid,'null,%d\n',i);
end
% fclose(fid);
