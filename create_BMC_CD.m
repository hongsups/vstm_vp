function [ML_estimate_mat, mean_L, sem_L, L_mat] = create_BMC_CD(expid,max_idx)

% clear all
% expid = 1;

file_path = 'CD/summary/';
if expid == 1       % orientation
    exp_name = 'orientation';
    nsbj = 10;
elseif expid ==2    % color
    exp_name = 'color';
    nsbj = 7;    
end

model_names = {'IL', 'SA', 'EP', 'VP'};
IL_models = {'F', 'U', 'P'};
SA_models = {'F', 'U', 'P'};
EP_models = {'X', 'F', 'U', 'P'};
VP_models = {'X', 'F', 'U', 'P'};
% n_param = [3 3 4 5];                % # of parameters per model (IL: 3, SA: 3, EP: 4, VP: 5)

% parameter ranges
BT_path = 'CD/bigtable/';
param_name = {'{\epsilon}', 'K', '{\lambda}', [], [] ; 'J1', 'prior', 'K', [], []; 'J', '{\alpha}', 'prior', 'K', []; 'J_{bar}', '{\tau}', '{\alpha}', 'prior', 'K'};
% SA: J1_vec, p_change_vec, capacity_idx
files = dir([BT_path 'CD_T_SA_1.mat']);
load([BT_path files.name],'J1_vec','p_change_vec','capacity_SA','Nvec');
% EP: J_EP_vec, power_vec, p_change_vec, capacity_vec
files = dir([BT_path 'CD_T_EP_1.mat']);
load([BT_path files.name],'J_EP_vec','power_vec');
% VP: Jbar_vec, tau_vec, power_vec, p_change_vec, capacity_vec;
files = dir([BT_path 'CD_T_VP_1_1.mat']);
load([BT_path files.name],'Jbar_vec','tau_vec','nsteps');
% IL: p_change_vec, capacity_vec, lambda
capacity_PK = linspace(1/3,10,nsteps);  % PK
capacity = 1:max(Nvec);                 % FK, UK
lambda = linspace(0, 1, nsteps);        
epsilon_vec = linspace(0, 1, nsteps); 

L_mat_sbj = [];
BIC_mat_sbj = [];

for i = 1:4
    model_id = model_names{i};
    capacity_vec = 0;
    if i == 1
        models = IL_models;
        param_range_mat = {epsilon_vec, capacity_vec, lambda};
    elseif i == 2
        models = SA_models;
        param_range_mat = {J1_vec, p_change_vec, capacity_vec};
    elseif i == 3
        models = EP_models;
        param_range_mat = {J_EP_vec, power_vec, p_change_vec, capacity_vec};
    elseif i == 4
        models = VP_models;
        param_range_mat = {Jbar_vec, tau_vec, power_vec, p_change_vec, capacity_vec};
    end
    
    for j = 1:length(models)
        model_name = [model_id '_' models{j}];
        files = dir([file_path 'LLH_' model_name '_' num2str(expid) '*.mat']);

        ML_estimate_mat_sbj = [];
        for k = 1:length(files)
            load([file_path files(k).name],'L','BIC','ML_estimate');
            
            % individual plot
            L_mat_sbj = [L_mat_sbj; L];
            BIC_mat_sbj = [BIC_mat_sbj; BIC];
            
            % ML estimates of a model from all subjects
            ML_estimate_mat_sbj = [ML_estimate_mat_sbj; ML_estimate];
        end
        
        ML_estimate_mat{i,j} = ML_estimate_mat_sbj;

        % boxplot for ML estimates
%         figure
%         for ii = 1:size(ML_estimate_mat_sbj,2)
%             if ii == (i + 1)
%                 if i > 2 && j < 4
%                     param_range_mat{ii} = capacity;
%                 elseif i == 1 && j < 3
%                     param_range_mat{ii} = capacity;
%                 elseif i == 2 && j < 3
%                     param_range_mat{ii} = capacity_SA;
%                 else
%                     param_range_mat{ii} = capacity_PK;
%                 end
%             end
%             subplot(1,size(ML_estimate_mat_sbj,2),ii);
%             boxplot(ML_estimate_mat_sbj(:,ii));title(model_name);xlabel([param_name{i,ii} ':' num2str(mean(ML_estimate_mat_sbj(:,ii)),3) '+/-' num2str(std((ML_estimate_mat_sbj(:,ii)))/sqrt(length((ML_estimate_mat_sbj(:,ii)))),3)]);
%             ylim([min(param_range_mat{ii}) max(param_range_mat{ii})])
%         end
    end
end

% compare mine and Shaiyan's
if expid == 1
    load all_max_params_orientation.mat
    ML_estimate_SK = all_max_params_orientation;
else
    load all_max_params_color.mat
    ML_estimate_SK = all_max_params_color;
end

% comparing SK's data with mine
for i = 1:4
    if i == 1       % IL (K, epsilon, lambda)
        figure;
        set(gca,'FontSize',11,'FontName','Arial');fac = 1;
        set(gcf,'Position',get(gcf,'Position').*[.1 .1 3 .8]*fac);
        set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);

        subplot(1,4,1); % K
        [corr_IL_K pValue_IL_K] = corrcoef(ML_estimate_mat{1,1}(:,2),ML_estimate_SK{1,1});
        plot(ML_estimate_SK{1,1}, ML_estimate_mat{1,1}(:,2)','o');title('K (IP-F)');
        hold on; plot(1:8,1:8,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_IL_K(2),3) ', p= ' num2str(pValue_IL_K(2),3)]);
        xlim([1 8]);ylim([1 8]);

        subplot(1,4,2); % epsilon
        [corr_IL_e pValue_IL_e] = corrcoef(ML_estimate_mat{1,1}(:,1),ML_estimate_SK{1,2});
        plot(ML_estimate_SK{1,2}, ML_estimate_mat{1,1}(:,1)','o');title('{\epsilon} (IP-F)');
        hold on; plot(0:.1:1,0:.1:1,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_IL_e(2),3) ', p= ' num2str(pValue_IL_e(2),3)]);
        xlim([0 1]);ylim([0 1]);

        subplot(1,4,3); % lambda
        [corr_IL_l pValue_IL_l] = corrcoef(ML_estimate_mat{1,1}(:,3),ML_estimate_SK{1,3});
        plot(ML_estimate_SK{1,3}, ML_estimate_mat{1,1}(:,3)','o');title('{\lambda} (IP-F)');
        hold on; plot(0:.1:1,0:.1:1,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_IL_l(2),3) ', p= ' num2str(pValue_IL_l(2),3)]);
        xlim([0 1]);ylim([0 1]);
        
        subplot(1,4,4); % errors
        boxplot([ML_estimate_SK{1,1}-ML_estimate_mat{1,1}(:,2)'; ML_estimate_SK{1,2}-ML_estimate_mat{1,1}(:,1)'; ML_estimate_SK{1,3}-ML_estimate_mat{1,1}(:,3)']');
        set(gca,'Xtick',1:3,'XtickLabel',{'K','e','g'});        
        
    elseif i == 2   % SA (J, K, p_change)
        figure;
        set(gca,'FontSize',11,'FontName','Arial');fac = 1;
        set(gcf,'Position',get(gcf,'Position').*[.1 .1 3 .8]*fac);
        set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);

        subplot(1,4,1); % J
        [corr_SA_J pValue_SA_J] = corrcoef(ML_estimate_mat{2,1}(:,1),ML_estimate_SK{2,1});
        plot(ML_estimate_SK{2,1}, ML_estimate_mat{2,1}(:,1)','o');title('J (SA-F)');
        hold on; plot(1:40,1:40,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_SA_J(2),3) ', p= ' num2str(pValue_SA_J(2),3)]);
        xlim([1 40]);ylim([1 40]);

        subplot(1,4,2); % K
        [corr_SA_K pValue_SA_K] = corrcoef(ML_estimate_mat{2,1}(:,3),ML_estimate_SK{2,2});
        plot(ML_estimate_SK{2,2}, ML_estimate_mat{2,1}(:,3)','o');title('K (SA-F)');
        hold on; plot(1:.1:25,1:.1:25,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_SA_K(2),3) ', p= ' num2str(pValue_SA_K(2),3)]);
        xlim([1 25]);ylim([1 25]);

        subplot(1,4,3); % p_change
        [corr_SA_pc pValue_SA_pc] = corrcoef(ML_estimate_mat{2,1}(:,2),ML_estimate_SK{2,3});
        plot(ML_estimate_SK{2,3}, ML_estimate_mat{2,1}(:,2)','o');title('P_{change} (SA-F)');
        hold on; plot(0:.1:1,0:.1:1,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_SA_pc(2),3) ', p= ' num2str(pValue_SA_pc(2),3)]);
        xlim([0 1]);ylim([0 1]);

        subplot(1,4,4); % errors
        boxplot([ML_estimate_SK{2,1}-ML_estimate_mat{2,1}(:,1)'; ML_estimate_SK{2,2}-ML_estimate_mat{2,1}(:,3)'; ML_estimate_SK{2,3}-ML_estimate_mat{2,1}(:,2)']');
        set(gca,'Xtick',1:3,'XtickLabel',{'J','K','P_change'});        
        
    elseif i == 3   % EP (J, alpha, p_change)
        figure;
        set(gca,'FontSize',11,'FontName','Arial');fac = 1;
        set(gcf,'Position',get(gcf,'Position').*[.1 .1 3 .8]*fac);
        set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);

        subplot(1,4,1); % J
        [corr_EP_J pValue_EP_J] = corrcoef(ML_estimate_mat{3,1}(:,1),ML_estimate_SK{4,1});
        plot(ML_estimate_SK{4,1}, ML_estimate_mat{3,1}(:,1)','o');title('J (EP-X)');
        hold on; plot(1:200,1:200,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_EP_J(2),3) ', p= ' num2str(pValue_EP_J(2),3)]);
        xlim([1 200]);ylim([1 200]);

        subplot(1,4,2); % alpha
        [corr_EP_a pValue_EP_a] = corrcoef(ML_estimate_mat{3,1}(:,3),-ML_estimate_SK{4,2});
        plot(-ML_estimate_SK{4,2}, ML_estimate_mat{3,1}(:,2)','o');title('{\alpha} (EP-X)');
        hold on; plot(0:.1:2,0:.1:2,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_EP_a(2),3) ', p= ' num2str(pValue_EP_a(2),3)]);
        xlim([0 2]);ylim([0 2]);

        subplot(1,4,3); % p_change
        [corr_EP_pc pValue_EP_pc] = corrcoef(ML_estimate_mat{3,1}(:,2),ML_estimate_SK{4,3});
        plot(ML_estimate_SK{4,3}, ML_estimate_mat{3,1}(:,3)','o');title('P_{change} (EP-X)');
        hold on; plot(0:.1:1,0:.1:1,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_EP_pc(2),3) ', p= ' num2str(pValue_EP_pc(2),3)]);
        xlim([0 1]);ylim([0 1]);

        subplot(1,4,4); % errors
        boxplot([ML_estimate_SK{4,1}-ML_estimate_mat{3,1}(:,1)'; ML_estimate_SK{4,2}-ML_estimate_mat{3,1}(:,2)'; ML_estimate_SK{4,3}-ML_estimate_mat{3,1}(:,3)']');
        set(gca,'Xtick',1:3,'XtickLabel',{'J','alpha','P_change'});        
        
    elseif i == 4   % VP (J, tau, alpha, p_change)
        figure;
        set(gca,'FontSize',11,'FontName','Arial');fac = 1;
        set(gcf,'Position',get(gcf,'Position').*[.1 .1 3 .8]*fac);
        set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);

        subplot(1,5,1); % J
        [corr_VP_J pValue_VP_J] = corrcoef(ML_estimate_mat{4,1}(:,1),ML_estimate_SK{5,1});
        plot(ML_estimate_SK{5,1}, ML_estimate_mat{4,1}(:,1)','o');title('J (VP-X)');
        hold on; plot(1:300,1:300,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_VP_J(2),3) ', p= ' num2str(pValue_VP_J(2),3)]);
        xlim([1 300]);ylim([1 300]);

        subplot(1,5,2); % tau
        [corr_VP_t pValue_VP_t] = corrcoef(ML_estimate_mat{4,1}(:,2),ML_estimate_SK{5,2});
        plot(ML_estimate_SK{5,2}, ML_estimate_mat{4,1}(:,2)','o');title('{\tau} (VP-X)');
        hold on; plot(0:.1:300,0:.1:300,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_VP_t(2),3) ', p= ' num2str(pValue_VP_t(2),3)]);
        xlim([0 300]);ylim([0 300]);        
        
        subplot(1,5,3); % alpha
        [corr_VP_a pValue_VP_a] = corrcoef(ML_estimate_mat{4,1}(:,3),-ML_estimate_SK{5,3});
        plot(-ML_estimate_SK{5,3}, ML_estimate_mat{4,1}(:,3)','o');title('{\alpha} (VP-X)');
        hold on; plot(0:.1:2,0:.1:2,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_VP_a(2),3) ', p= ' num2str(pValue_VP_a(2),3)]);
        xlim([0 2]);ylim([0 2]);

        subplot(1,5,4); % p_change
        [corr_VP_pc pValue_VP_pc] = corrcoef(ML_estimate_mat{4,1}(:,4),ML_estimate_SK{5,4});
        plot(ML_estimate_SK{5,4}, ML_estimate_mat{4,1}(:,4)','o');title('P_{change} (VP-X)');
        hold on; plot(0:.1:1,0:.1:1,':');ylabel('HS');
        xlabel(['corr= ' num2str(corr_VP_pc(2),3) ', p= ' num2str(pValue_VP_pc(2),3)]);
        xlim([0 1]);ylim([0 1]);
        
        subplot(1,5,5); % errors
        boxplot([ML_estimate_SK{5,1}-ML_estimate_mat{4,1}(:,1)'; ML_estimate_SK{5,2}-ML_estimate_mat{4,1}(:,2)'; ML_estimate_SK{5,3}-ML_estimate_mat{4,1}(:,3)'; ML_estimate_SK{5,4}-ML_estimate_mat{4,1}(:,4)']');
        set(gca,'Xtick',1:3,'XtickLabel',{'J','tau','alpha','P_change'});        
       
    end
end

% L_mat_sbj has every subject LLH value -> need to be re-shaped
L_mat = reshape(L_mat_sbj,nsbj,14);     % 14 = model number
BIC_mat = reshape(BIC_mat_sbj,nsbj,14);     % 14 = model number

load modelnames_14_name_final.mat;
m = [sbjname_IL; sbjname_SA; sbjname_EP; sbjname_VP];
% max_idx = 14;
max_model_name = m{max_idx};
L_relative = L_mat - repmat(L_mat(:, max_idx),1,length(m));
mean_L = mean(L_relative);
sem_L = std(L_relative)/sqrt(nsbj);

% models x subject (3 or 4 x nsbj): encoding precision
L_mat_IL = L_mat(:,1:3)';
L_mat_SA = L_mat(:,4:6)';
L_mat_EP = L_mat(:,7:10)';
L_mat_VP = L_mat(:,11:14)';

% reference point
L_mat_ref = L_mat(:,max_idx)';

L_mat_IL_subtracted = L_mat_IL - repmat(L_mat_ref,3,1);
L_mat_SA_subtracted = L_mat_SA - repmat(L_mat_ref,3,1);
L_mat_EP_subtracted = L_mat_EP - repmat(L_mat_ref,4,1);
L_mat_VP_subtracted = L_mat_VP - repmat(L_mat_ref,4,1);

% L_mat_ref = max(L_mat(:,max_idx));

% L_mat_IL_subtracted = L_mat_IL - L_mat_ref;
% L_mat_SA_subtracted = L_mat_SA - L_mat_ref;
% L_mat_EP_subtracted = L_mat_EP - L_mat_ref;
% L_mat_VP_subtracted = L_mat_VP - L_mat_ref;

% 1. create nbsj * 4 
L_IL = log(mean(exp(L_mat_IL_subtracted)));
L_SA = log(mean(exp(L_mat_SA_subtracted)));
L_EP = log(mean(exp(L_mat_EP_subtracted)));
L_VP = log(mean(exp(L_mat_VP_subtracted)));

L_IL_mean_over_sbj = mean(L_IL - L_VP);
L_SA_mean_over_sbj = mean(L_SA - L_VP);
L_EP_mean_over_sbj = mean(L_EP - L_VP);
L_VP_mean_over_sbj = mean(L_VP - L_VP);

sem_IL = std(L_IL - L_VP)/sqrt(nsbj);
sem_SA = std(L_SA - L_VP)/sqrt(nsbj);
sem_EP = std(L_EP - L_VP)/sqrt(nsbj);
sem_VP = std(L_VP - L_VP)/sqrt(nsbj);

% L_mat_IL_subtracted = L_mat_IL - repmat(max(L_mat_IL),3,1);
% L_mat_SA_subtracted = L_mat_SA - repmat(max(L_mat_SA),3,1);
% L_mat_EP_subtracted = L_mat_EP - repmat(max(L_mat_EP),4,1);
% L_mat_VP_subtracted = L_mat_VP - repmat(max(L_mat_VP),4,1);
% 
% L_mat_IL_mean = min(L_mat_IL) - log(nsbj) + log(sum(exp(L_mat_IL_subtracted)));
% L_mat_SA_mean = min(L_mat_SA) - log(nsbj) + log(sum(exp(L_mat_SA_subtracted)));
% L_mat_EP_mean = min(L_mat_EP) - log(nsbj) + log(sum(exp(L_mat_EP_subtracted)));
% L_mat_VP_mean = min(L_mat_VP) - log(nsbj) + log(sum(exp(L_mat_VP_subtracted)));
% 
% L_IL_mean_over_sbj = mean(L_mat_IL_mean-L_mat_VP_mean);
% L_SA_mean_over_sbj = mean(L_mat_SA_mean-L_mat_VP_mean);
% L_EP_mean_over_sbj = mean(L_mat_EP_mean-L_mat_VP_mean);
% L_VP_mean_over_sbj = mean(L_mat_VP_mean-L_mat_VP_mean);
% 
% sem_IL = std(L_mat_IL_mean-L_mat_VP_mean)/sqrt(nsbj);
% sem_SA = std(L_mat_SA_mean-L_mat_VP_mean)/sqrt(nsbj);
% sem_EP = std(L_mat_EP_mean-L_mat_VP_mean)/sqrt(nsbj);
% sem_VP = std(L_mat_VP_mean-L_mat_VP_mean)/sqrt(nsbj);

L_group_mean = [L_IL_mean_over_sbj L_SA_mean_over_sbj L_EP_mean_over_sbj L_VP_mean_over_sbj];
sem_group_mean = [sem_IL sem_SA sem_EP sem_VP];

% group mean: # of encoded items
L_mat_A = L_mat(:,[7 11])';
L_mat_F = L_mat(:,[1 4 8 12])';
L_mat_P = L_mat(:,[3 6 10 14])';
L_mat_U = L_mat(:,[2 5 9 13])';

L_mat_A_mean = min(L_mat_A) - log(nsbj) + log(sum(exp(L_mat_A - repmat(max(L_mat_A),2,1))));
L_mat_F_mean = min(L_mat_F) - log(nsbj) + log(sum(exp(L_mat_F - repmat(max(L_mat_F),4,1))));
L_mat_P_mean = min(L_mat_P) - log(nsbj) + log(sum(exp(L_mat_P - repmat(max(L_mat_P),4,1))));
L_mat_U_mean = min(L_mat_U) - log(nsbj) + log(sum(exp(L_mat_U - repmat(max(L_mat_U),4,1))));

L_A_mean_over_sbj = mean(L_mat_A_mean-L_mat_P_mean);
L_F_mean_over_sbj = mean(L_mat_F_mean-L_mat_P_mean);
L_P_mean_over_sbj = mean(L_mat_P_mean-L_mat_P_mean);
L_U_mean_over_sbj = mean(L_mat_U_mean-L_mat_P_mean);

sem_A = std(L_mat_A_mean-L_mat_P_mean)/sqrt(nsbj);
sem_F = std(L_mat_F_mean-L_mat_P_mean)/sqrt(nsbj);
sem_P = std(L_mat_P_mean-L_mat_P_mean)/sqrt(nsbj);
sem_U = std(L_mat_U_mean-L_mat_P_mean)/sqrt(nsbj);

L_group_mean_K = [L_A_mean_over_sbj L_F_mean_over_sbj L_U_mean_over_sbj L_P_mean_over_sbj];
sem_group_mean_K = [sem_A sem_F sem_U sem_P];

L_sorted = sort(mean_L,'ascend');
for i = 1:length(L_sorted)
    m_sorted{i} = m{L_sorted(i)==mean_L};
    sem_L_sorted(i) = sem_L(L_sorted(i)==mean_L);
end
close all

% average BMC
figure;
set(gca,'FontSize',11,'FontName','Arial');
fac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 2 2]*fac);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);
barh(L_sorted,'FaceColor','y');
xlabel(['Model LLH Relative to ' char(max_model_name) ' Model']);xlim([-160 10]);ylim([0 15]);ylabel('Models');title(['Bayseian model comparison on average subject data (' exp_name ')']);
hold all
errorbarh(1:length(m),L_sorted,sem_L_sorted');
set(gca,'YAxisLocation','right')
set(gca,'Ytick',1:length(m),'YtickLabel',m_sorted);
saveas(gcf,['BMC_CD_' exp_name '.emf'],'emf');

%% old one (encoding precision)
figure;
set(gca,'FontSize',11,'FontName','Arial');
fac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 2 2]*fac);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);
bar(L_group_mean,'FaceColor','y');
% xlabel(['Model LLH Relative to ' char(max_model_name) ' Model']);
% xlim([-160 10]);ylim([0 15]);
ylabel('Model likelihood');
title(['Bayseian model comparison on average subject data (encoding precision, ' exp_name ')']);
hold all
h = errorbar(1:length(L_group_mean),L_group_mean,sem_group_mean);
set(h,'linestyle','none');
set(gca,'Xtick',1:length(L_group_mean));
set(gca,'XTickLabel',{'IP','SA','EP','VP'});
% set(gca,'Ytick',-10:2:0);
saveas(gcf,['BMC_CD_group_precision_' exp_name '.emf'],'emf');

figure;
set(gca,'FontSize',11,'FontName','Arial');
fac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 2 2]*fac);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);
bar(L_group_mean_K,'FaceColor','y');
% xlabel(['Model LLH Relative to ' char(max_model_name) ' Model']);
% xlim([-160 10]);ylim([0 15]);
ylabel('Model likelihood');
title(['Bayseian model comparison on average subject data (# of items, ' exp_name ')']);
hold all
h = errorbar(1:length(L_group_mean_K),L_group_mean_K,sem_group_mean_K);
set(h,'linestyle','none');
set(gca,'Xtick',1:length(L_group_mean_K));
set(gca,'XTickLabel',{'All','Fixed','Uniform','Poisson'});
% set(gca,'Ytick',-10:2:0);
saveas(gcf,['BMC_CD_group_K_' exp_name '.emf'],'emf');

% % individual BMC
% figure;
% set(gca,'FontSize',11,'FontName','Arial');
% fac = 1;
% set(gcf,'Position',get(gcf,'Position').*[.1 .1 2 2]*fac);
% set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);
% bar(L_relative);

% xlabel(['Model LLH Relative to ' char(max_model_name) ' Model']);xlim([-160 10]);ylim([0 15]);ylabel('Models');title(['Bayseian model comparison on average subject data (' exp_name ')']);

function errorbarh(Y,X,stderr)
for ii=1:length(Y)
    plot([X(ii)-stderr(ii) X(ii)+stderr(ii)],[Y(ii) Y(ii)],'k');
end
