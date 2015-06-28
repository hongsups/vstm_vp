function [L, ML_estimate, BIC] = create_LLH_IL_CD(capacity_idx,expid,subjid)

if expid < 3
    file_path = 'CD/summary/';
    [~, model_name] = create_modelname([],[],[],[],capacity_idx);
    sbj_name = [num2str(expid) 'S' num2str(subjid)];
    fname = [file_path 'LLH_IL_' model_name '_' sbj_name '.mat'];
    fname2 = [file_path 'perf_IL_' model_name '_' sbj_name '.mat'];
else
    file_path = 'CD/faketest/';
    model_name = create_modelname_temp(1,capacity_idx);
    fname = [file_path 'LLH_' model_name '_' subjid '.mat'];
    fname2 = [file_path 'perf_' model_name '_' subjid '.mat'];
end

% if ~exist(fname,'file')

    Nvec = 1:8;                             % set sizes
    if expid == 1
        Nvec_exp = 2:2:8;                       % set sizes used in the experiments
    else
        Nvec_exp = [1 2 4 8];
    end
    sampleSize_n = 500;                     % sample size
    nsteps = 20;                            % # of steps in a parameter space

    capacity_PK = linspace(1/3,10,nsteps);  % PK
    capacity = 1:max(Nvec);                 % FK, UK
    % When no change occurs among the memorized items, the observer responds ?change? with a guessing probability g.
    lambda = linspace(0, 1, nsteps);        
    % When a change occurs among the memorized items, the observer responds ?change? with probability 1-e.
    epsilon_vec = linspace(0, 1, nsteps); 

    if capacity_idx == 4
        capacity_vec = capacity_PK;
    else
        capacity_vec = capacity;
    end

    if capacity_idx == 1
        error('Capacity cannot be inifinite.');
    end 

    % create look-up table
    p_res_change = zeros(2,length(epsilon_vec),length(capacity),length(lambda),length(Nvec));
    for C = 1:2
        for ie = 1:length(epsilon_vec)
            epsilon = epsilon_vec(ie);
            for ic = 1:length(capacity)
                K = capacity(ic);
                for il = 1:length(lambda)
                    g = lambda(il);
                    for n = 1:length(Nvec)
                        N = Nvec(n);
                        if C == 1       % change
                            if N>K
                                p_res_change(C, ie, ic, il, n) = K/N*(1-epsilon)+(1-K/N)*g;
                            else
                                p_res_change(C, ie, ic, il, n) = 1-epsilon;
                            end
                        else            % no change
                            p_res_change(C, ie, ic, il, n) = g;
                        end                        
                    end
                end
            end
        end
    end

    % load subject data
    data = loaddata_CD(subjid,expid);
    response_idx_vec = data(:,2);
    N_idx_vec = data(:,5);
%     if expid < 3
%         change_idx_vec = ceil((data(:,1))/max(data(:,1)));
%     else
%         change_idx_vec = data(:,1);
%     end
    if expid < 3
        change_idx_vec = ceil((data(:,1))/max(data(:,1)));
        mag_vec = data(:,6)/pi*180.*change_idx_vec;     % [0,180] (You need to multiply change_idx)    
    else
        change_idx_vec = data(:,1);
        mag_vec = abs(data(:,6));    % For fake data, you don't have to multiply change_idx_vec. It's already take care of in fake data generatio code.
    end
    % categorizing trials by N
    hit_trials = zeros(length(Nvec_exp),1);
    false_alarm_trials = zeros(length(Nvec_exp),1);
    miss_trials = zeros(length(Nvec_exp),1);
    correct_rejection_trials = zeros(length(Nvec_exp),1);
    for i = 1:length(Nvec_exp)
        hit_trials(i) = sum((response_idx_vec == 1) & (change_idx_vec == 1) & (N_idx_vec == Nvec_exp(i)));              % HIT
        false_alarm_trials(i) = sum((response_idx_vec == 1) & (change_idx_vec == 0) & (N_idx_vec == Nvec_exp(i)));      % FALSE ALARM
        miss_trials(i) = sum((response_idx_vec == 0) & (change_idx_vec == 1) & (N_idx_vec == Nvec_exp(i)));             % MISS
        correct_rejection_trials(i) = sum((response_idx_vec == 0) & (change_idx_vec == 0) & (N_idx_vec == Nvec_exp(i)));% CORRECT REJECTION
    end

    % precompute weights for poisson-K model
    maxKinPoisson = 25;   % we cut off the poisson distribution here, because hardly any prob mass left
    poiss_w = zeros(length(capacity_PK),9);
    for ii=1:length(capacity_PK)
        y = poisspdf(0:maxKinPoisson,capacity_PK(ii));
        poiss_w(ii,1:8) = y(1:8);
        poiss_w(ii,9) = sum(y(9:end));
        poiss_w(ii,:) = poiss_w(ii,:)/sum(poiss_w(ii,:));
    end

    % compute likelihood
    p_L = zeros(length(epsilon_vec),length(capacity_vec),length(lambda));
    if capacity_idx == 2                    % fixed K
        for ie = 1:length(epsilon_vec)
            for ic = 1:length(capacity)
                for il = 1:length(lambda)
                    p_response = zeros(2,length(Nvec_exp));
                    for i = 1:length(Nvec_exp)
                        p_response(:,i) = p_res_change(:,ie,ic,il,Nvec_exp(i));
%                         p_response(:,i) = p_res_change(:,ie,ic,il,i*2);
                    end
                    p_no_response = 1 - p_response;

                    p_response(p_response==0) = 1/(sampleSize_n*10);
                    p_no_response(p_no_response==0) = 1/(sampleSize_n*10);

                    p_log_like_hit = hit_trials'.*log(p_response(1,:));
                    p_log_like_false_alarm = false_alarm_trials'.*log(p_response(2,:));
                    p_log_like_miss = miss_trials'.*log(p_no_response(1,:));
                    p_log_like_correct_rejection = correct_rejection_trials'.*log(p_no_response(2,:));

                    p_log_like_mat = p_log_like_hit + p_log_like_false_alarm + p_log_like_miss + p_log_like_correct_rejection; 
                    p_L(ie,ic,il) = sum(p_log_like_mat(:));

                    
                end
            end
        end
        
    elseif capacity_idx == 3                % uniform K
        for ie = 1:length(epsilon_vec)
            for ic = 1:length(capacity)
                Kmax = capacity(ic);
                for il = 1:length(lambda)
                    p_response = zeros(2,length(Nvec_exp));
                    for i = 1:length(Nvec_exp)
                        for kk=1:Kmax
                            p_response(:,i) = p_response(:,i) + 1/Kmax*p_res_change(:,ie,kk,il,Nvec_exp(i));
%                             p_response(:,i) = p_response(:,i) + 1/Kmax*p_res_change(:,ie,kk,il,i*2);
                        end
                    end
                    p_no_response = 1-p_response;

                    p_response(p_response==0) = 1/(sampleSize_n*10);
                    p_no_response(p_no_response==0) = 1/(sampleSize_n*10);

                    p_log_like_hit = hit_trials'.*log(p_response(1,:));
                    p_log_like_false_alarm = false_alarm_trials'.*log(p_response(2,:));
                    p_log_like_miss = miss_trials'.*log(p_no_response(1,:));
                    p_log_like_correct_rejection = correct_rejection_trials'.*log(p_no_response(2,:));

                    p_log_like_mat = p_log_like_hit + p_log_like_false_alarm + p_log_like_miss + p_log_like_correct_rejection; 
                    p_L(ie,ic,il) = sum(p_log_like_mat(:));
                end
            end
        end
    elseif capacity_idx == 4                  % poisson K
        for ie = 1:length(epsilon_vec)
            for ic = 1:length(capacity_PK)
                for il = 1:length(lambda)
                    p_response = zeros(2,length(Nvec_exp));
                    p_response_for_zero_K = 1/2;
                    for i = 1:length(Nvec_exp)
                        % case K=0
                            p_response(:,i) = p_response(:,i) + poiss_w(ic,1)*p_response_for_zero_K;
                        % case K>0
                        for kk=1:8
%                             p_response(:,i) = p_response(:,i) + poiss_w(ic,kk+1)*p_res_change(:,ie,kk,il,i*2);
                            p_response(:,i) = p_response(:,i) + poiss_w(ic,kk+1)*p_res_change(:,ie,kk,il,Nvec_exp(i));
                        end                                    
%                         for kk=1:8
%                             p_response(:,i) = p_response(:,i) + poiss_w(ic,kk)*p_res_change(:,ie,kk,il,i*2);
%                         end
                    end
                    p_response = min(p_response,1);
                    p_no_response = 1-p_response;

                    p_response(p_response==0) = 1/(sampleSize_n*10);
                    p_no_response(p_no_response==0) = 1/(sampleSize_n*10);

                    p_log_like_hit = hit_trials'.*log(p_response(1,:));
                    p_log_like_false_alarm = false_alarm_trials'.*log(p_response(2,:));
                    p_log_like_miss = miss_trials'.*log(p_no_response(1,:));
                    p_log_like_correct_rejection = correct_rejection_trials'.*log(p_no_response(2,:));

                    p_log_like_mat = p_log_like_hit + p_log_like_false_alarm + p_log_like_miss + p_log_like_correct_rejection; 
                    p_L(ie,ic,il) = sum(p_log_like_mat(:));
                end
            end
        end
    end

    % model likelihood
    p_L_max = max(p_L(:));
    expL = exp(p_L-p_L_max);            

    % old
    % L = log(sum(expL(:))) + p_L_max - log(length(epsilon_vec)*length(capacity_vec)*length(lambda));

    % new
    L = log(trapz(epsilon_vec,trapz(capacity_vec,trapz(lambda,expL,3),2),1)) + p_L_max - log(range(epsilon_vec)*range(capacity_vec)*range(lambda));

    % ML estimate & BIC
    [M, I] = max(p_L(:));     
    [idx_e, idx_c, idx_l] = ind2sub(size(p_L),I);
    ML_estimate = [epsilon_vec(idx_e) capacity_vec(idx_c) lambda(idx_l)];
    BIC = max(p_L(:)) - 3/2*log(length(data));

    epsilon_est = ML_estimate(1);
    K_est = ML_estimate(2);
    g_est = ML_estimate(3);

    % model fitting
    nBins = 9;
    bin = linspace(0,90,nBins+1);
    bin = bin(1:end-1)+diff(bin(1:2))/2;
    binsize = diff(bin(1:2));

    % hit, false alarm
    hit_rate = zeros(1,length(Nvec_exp));
    false_alarm_rate = zeros(1,length(Nvec_exp));
    prop_magN = zeros(length(Nvec_exp),length(bin));
    m_hit_rate = zeros(1,length(Nvec_exp));
    m_false_alarm_rate = zeros(1,length(Nvec_exp));
    for i = 1:length(Nvec_exp)
        N = Nvec_exp(i);
        hit_rate(i) = sum((response_idx_vec == 1) & (change_idx_vec == 1) & (N_idx_vec == N))/sum(change_idx_vec == 1 & N_idx_vec == N);
        false_alarm_rate(i) = sum((response_idx_vec == 1) & (change_idx_vec == 0) & (N_idx_vec == N))/sum(change_idx_vec == 0 & N_idx_vec == N);
        for jj = 1:length(bin)
            ind_sbj = (abs(mag_vec) >= bin(jj)-binsize/2) & (abs(mag_vec) < bin(jj)+binsize/2) & N_idx_vec == N & response_idx_vec==1 & change_idx_vec == 1;
            prop_magN(i,jj) = sum(ind_sbj)/sum(N_idx_vec == N & abs(mag_vec) >= bin(jj)-binsize/2 & abs(mag_vec) < bin(jj)+binsize/2 & change_idx_vec == 1);
        end    

        if N > K_est
            m_hit_rate(i) = K_est/N*(1-epsilon_est) + (1-K_est/N)*g_est;
        else
            m_hit_rate(i) = 1 - epsilon_est;
        end
        m_false_alarm_rate(i) = g_est;
    end
    m_prop_magN = repmat(m_hit_rate',1,length(bin));

    figure;
    set(gca,'FontSize',11,'FontName','Arial');
    fac = 1;
    set(gcf,'Position',get(gcf,'Position').*[.1 .1 2 2]*fac);
    set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);
    % performance by set size
    subplot(1,2,1); 
    plot(Nvec_exp, hit_rate,'ob');
    hold on;
    plot(Nvec_exp, false_alarm_rate,'or');
    hold on;
    plot(Nvec_exp, m_hit_rate,':b');
    hold on;
    plot(Nvec_exp, m_false_alarm_rate,':r');
    % plot(Nvec_exp, m_prop_response_N,':'); 
    ylabel('Proportion reports "change"');axis([2 8 0 1]);
    legend('hit rate','false alarm rate',1);

    % performance by change magnitude: add 0 on x-axis for FA rate
    bin_new = [0 bin];
    prop_magN_new = [false_alarm_rate' prop_magN];
    m_prop_magN_new = [m_false_alarm_rate' m_prop_magN];

    subplot(1,2,2); 
    plot(bin_new, prop_magN_new, 'o');
    hold on
    plot(bin_new, m_prop_magN_new, ':'); ylabel('Proportion correct');
    axis([0 90 0 1]);legend(strcat('N= ',int2str(Nvec_exp')), 'Location','Best');

    save(fname,'L','ML_estimate','BIC');       
    save(fname2,'hit_rate','false_alarm_rate','prop_magN','m_hit_rate','m_false_alarm_rate','m_prop_magN','bin','bin_new','prop_magN_new','m_prop_magN_new');

% else
%     fprintf('File already exists - skipping\n');
% end
