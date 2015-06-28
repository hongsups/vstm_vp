function create_modelfit_CD(expid,modelidx)

% model fit for subject data
file_path = 'CD/summary/';
if expid == 1       % orientation
    exp_name = 'ori';
    nsbj = 10;
elseif expid ==2    % color
    exp_name = 'col';
    nsbj = 7;    
end

if expid == 1
    Nvec = 2:2:8;                       % set sizes used in the experiments
else
    Nvec = [1 2 4 8];
end
nBins = 9;
bin = linspace(0,90,nBins+1);
bin = bin(1:end-1)+diff(bin(1:2))/2;
binsize = diff(bin(1:2));
colorvec = get(gca, 'ColorOrder');
colorvec = min(colorvec+.65,1);
close;

load modelnames_14_name_final.mat;
capacity_idx_vec = 'AFUP';
model_names = {'IL', 'SA', 'EP', 'VP'};
ntrials = 5000;
%for i = [7 11] % 1:length(sbjname) % original

i = modelidx;
model_name = sbjname{i};
capacity_idx = find(capacity_idx_vec == model_name(3));    

% we need to make a fake data...
if i < 4
    files = dir([file_path 'perf_IL_' model_name(3) '_' num2str(expid) 'S*.mat']);

    hit_rate_vec = [];
    false_alarm_rate_vec = [];
    prop_magN_vec = [];
    m_hit_rate_vec = [];
    m_false_alarm_rate_vec = [];
    m_prop_magN_vec = [];
    for j = 1:length(files)
        load([file_path files(j).name]);
        hit_rate_vec = [hit_rate_vec; hit_rate];
        false_alarm_rate_vec = [false_alarm_rate_vec; false_alarm_rate];
        prop_magN_vec = [prop_magN_vec; prop_magN];
        m_hit_rate_vec = [m_hit_rate_vec; m_hit_rate];
        m_false_alarm_rate_vec = [m_false_alarm_rate_vec; m_false_alarm_rate];
        m_prop_magN_vec = [m_prop_magN_vec; m_prop_magN];
    end
else
    if sum(i == 4:6) == 1
        files = dir([file_path 'LLH_SA_' model_name(3) '_' num2str(expid) 'S*.mat']);

        hit_rate_vec = [];
        false_alarm_rate_vec = [];
        prop_magN_vec = [];
        m_hit_rate_vec = [];
        m_false_alarm_rate_vec = [];
        m_prop_magN_vec = [];
        for j = 1:length(files)
            load([file_path files(j).name]);
            Jbar = ML_estimate(1);
            p_change = ML_estimate(2);
            K = ML_estimate(3);
            data_model = create_fakedata_SA_CD(ntrials,capacity_idx,Jbar,p_change,K,expid);
            data_sbj = loaddata_CD(j, expid);
            [hit_rate false_alarm_rate prop_magN m_hit_rate m_false_alarm_rate m_prop_magN] = create_fit_CD(data_sbj,data_model,expid);
            hit_rate_vec = [hit_rate_vec; hit_rate];
            false_alarm_rate_vec = [false_alarm_rate_vec; false_alarm_rate];
            prop_magN_vec = [prop_magN_vec; prop_magN];
            m_hit_rate_vec = [m_hit_rate_vec; m_hit_rate];
            m_false_alarm_rate_vec = [m_false_alarm_rate_vec; m_false_alarm_rate];
            m_prop_magN_vec = [m_prop_magN_vec; m_prop_magN];
        end
    elseif sum(i == 7:10) == 1
        if strcmp(model_name(3),'A')            
            files = dir([file_path 'LLH_EP_X_' num2str(expid) 'S*.mat']);
        else
            files = dir([file_path 'LLH_EP_' model_name(3) '_' num2str(expid) 'S*.mat']);
        end

        hit_rate_vec = [];
        false_alarm_rate_vec = [];
        prop_magN_vec = [];
        m_hit_rate_vec = [];
        m_false_alarm_rate_vec = [];
        m_prop_magN_vec = [];
        for j = 1:length(files)
            load([file_path files(j).name]);

            Jbar = ML_estimate(1);
            power = ML_estimate(2);
            p_change = ML_estimate(3);
            K = ML_estimate(4);

            data_model = create_fakedata_EP_CD(ntrials,capacity_idx,Jbar,power,p_change,K,expid);
            data_sbj = loaddata_CD(j, expid);
            [hit_rate false_alarm_rate prop_magN m_hit_rate m_false_alarm_rate m_prop_magN] = create_fit_CD(data_sbj,data_model,expid);

            hit_rate_vec = [hit_rate_vec; hit_rate];
            false_alarm_rate_vec = [false_alarm_rate_vec; false_alarm_rate];
            prop_magN_vec = [prop_magN_vec; prop_magN];
            m_hit_rate_vec = [m_hit_rate_vec; m_hit_rate];
            m_false_alarm_rate_vec = [m_false_alarm_rate_vec; m_false_alarm_rate];
            m_prop_magN_vec = [m_prop_magN_vec; m_prop_magN];
        end

    elseif sum(i == 11:14) == 1        
        if strcmp(model_name(3),'A')            
            files = dir([file_path 'LLH_VP_X_' num2str(expid) 'S*.mat']);
        else
            files = dir([file_path 'LLH_VP_' model_name(3) '_' num2str(expid) 'S*.mat']);
        end

        hit_rate_vec = [];
        false_alarm_rate_vec = [];
        prop_magN_vec = [];
        m_hit_rate_vec = [];
        m_false_alarm_rate_vec = [];
        m_prop_magN_vec = [];
        for j = 1:length(files)
            load([file_path files(j).name]);

            Jbar = ML_estimate(1);
            t = ML_estimate(2);
            power = ML_estimate(3);
            p_change = ML_estimate(4);
            K = ML_estimate(5);

            data_model = create_fakedata_VP_CD(ntrials,capacity_idx,Jbar,t,power,p_change,K,expid);
            data_sbj = loaddata_CD(j, expid);
            [hit_rate false_alarm_rate prop_magN m_hit_rate m_false_alarm_rate m_prop_magN] = create_fit_CD(data_sbj,data_model,expid);

            hit_rate_vec = [hit_rate_vec; hit_rate];
            false_alarm_rate_vec = [false_alarm_rate_vec; false_alarm_rate];
            prop_magN_vec = [prop_magN_vec; prop_magN];
            m_hit_rate_vec = [m_hit_rate_vec; m_hit_rate];
            m_false_alarm_rate_vec = [m_false_alarm_rate_vec; m_false_alarm_rate];
            m_prop_magN_vec = [m_prop_magN_vec; m_prop_magN];
        end
    end
end

mean_hit_rate_vec = mean(hit_rate_vec);
mean_false_alarm_rate_vec = mean(false_alarm_rate_vec);
for n = 1:length(Nvec)
    mean_prop_magN_vec(n, :) = mean(prop_magN_vec(n:length(Nvec):end, :));
    mean_m_prop_magN_vec(n, :) = mean(m_prop_magN_vec(n:length(Nvec):end, :));
    sem_prop_magN_vec(n, :) = std(prop_magN_vec(n:length(Nvec):end, :))/sqrt(nsbj);
    sem_m_prop_magN_vec(n, :) = std(m_prop_magN_vec(n:length(Nvec):end, :))/sqrt(nsbj);
end
mean_m_hit_rate_vec = mean(m_hit_rate_vec);
mean_m_false_alarm_rate_vec = mean(m_false_alarm_rate_vec);

sem_hit_rate_vec = std(hit_rate_vec)/sqrt(nsbj);
sem_false_alarm_rate_vec = std(false_alarm_rate_vec)/sqrt(nsbj);
sem_m_hit_rate_vec = std(m_hit_rate_vec)/sqrt(nsbj);
sem_m_false_alarm_rate_vec = std(m_false_alarm_rate_vec)/sqrt(nsbj);

RMSE = sqrt(mean((prop_magN_vec(:)-m_prop_magN_vec(:)).^2));

%     % plotting
%     figure;        
%     set(gca,'FontSize',11,'FontName','Arial');
%     fac = .7;
%     set(gcf,'Position',get(gcf,'Position').*[.1 .1 1 1]*fac);
%     set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 1 1]*fac);

% performance by change magnitude: add 0 on x-axis for FA rate
bin_new = [0 bin];
if strcmp(model_name(1),'1')        % IL model
    mean_m_prop_magN_vec = repmat(mean_m_hit_rate_vec',1,length(bin));
    sem_m_prop_magN_vec = repmat(sem_m_hit_rate_vec',1,length(bin));
end
close all

prop_magN_new = [mean_false_alarm_rate_vec' mean_prop_magN_vec];
sem_prop_magN_new = [zeros(4,1) sem_prop_magN_vec];
m_prop_magN_new = [mean_m_false_alarm_rate_vec' mean_m_prop_magN_vec];
sem_m_prop_magN_new = [zeros(4,1) sem_m_prop_magN_vec];

% model fit
for pp = 1:length(Nvec)
    patch([bin_new, wrev(bin_new)],[m_prop_magN_new(pp, :) - sem_m_prop_magN_new(pp, :), wrev(m_prop_magN_new(pp, :) + sem_m_prop_magN_new(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;
% subject data
errorbar(repmat(bin_new, 4, 1)', prop_magN_new', sem_prop_magN_new', 'o');
% xlabel('Change Magnitude ({\Delta})');ylabel('P(C_{hat}=1)');
title([model_name '-' exp_name]);
axis([0 90 0 1]);
set(gca,'YTick',0:.2:1);set(gca,'XTick',0:30:90)
text(20,.1,['RMSE' '=' num2str(RMSE, 3)]);
% legend(strcat('N= ',int2str(Nvec')), 'Location','Best');

saveas(gcf,['fit_CD_' model_name '_' num2str(expid) '.emf'],'emf');

% save all
% qdas

% saveas(gcf,['modelfit_' model_id '_' exp 'CD.emf'],'emf');



