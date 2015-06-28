function [L, ML_estimate, BIC] = create_LLH_VP_CD(capacity_idx,expid,subjid)

% load parameters from the table file
BT_path = 'CD/bigtable/';
files = dir([BT_path 'CD_T_VP_1_1.mat']);
load([BT_path files.name],'nsteps','Jbar_vec','power_vec','tau_vec','p_change_vec','magstep','mag','sampleSize_n','Nvec');

if expid == 1
    Nvec_exp = 2:2:8;                       % set sizes used in the experiments
else
    Nvec_exp = [1 2 4 8];
end

capacity_PK = linspace(1/3,10,nsteps);  % PK
capacity = 1:max(Nvec);                 % FK, UK

if capacity_idx == 4
    capacity_vec = capacity_PK;
else
    capacity_vec = capacity;
end

% load subject data
data = loaddata_CD(subjid,expid);
if expid < 3
    change_idx_vec = ceil((data(:,1))/max(data(:,1)));
    mag_vec = data(:,6)/pi*180.*change_idx_vec;     % [0,180] (You need to multiply change_idx)    
    mag_vec(mag_vec>90)=180-mag_vec(mag_vec>90);    % ceil(mag_vec) -> min: 0, max: 90
else
    mag_vec = abs(data(:,6));   % For fake data, you don't have to multiply change_idx_vec. It's already take care of in fake data generatio code.
end
response_idx_vec = data(:,2);
N_idx_vec = data(:,5);

% remapping the delta vector to [0,90]
magidx = ceil(mag_vec/magstep)+1;
mag_s = 1:max(magidx);

% categorizing trials by N and magnitude
response_yes_trials = zeros(length(Nvec_exp),length(mag_s));
response_no_trials = zeros(length(Nvec_exp),length(mag_s));
for i = 1:length(Nvec_exp)
    N = Nvec_exp(i);
    for jj = 1:length(mag_s)        
        response_yes_trials(i,jj) = sum((response_idx_vec == 1) & (N_idx_vec == N) & (magidx == mag_s(jj)));     % Response: 'yes'
        response_no_trials(i,jj) = sum((response_idx_vec == 0) & (N_idx_vec == N) & (magidx == mag_s(jj)));      % Response: 'no'
    end
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

% p_res_change(ip, ipc, ic, i, j)

% compute likelihood
% initialize
if capacity_idx == 1
    p_L = zeros(length(Jbar_vec),length(tau_vec),length(power_vec),length(p_change_vec));
else
    p_L = zeros(length(Jbar_vec),length(tau_vec),length(power_vec),length(p_change_vec),length(capacity_vec));
end
for ij = 1:length(Jbar_vec)
    for it = 1:length(tau_vec)
        files = dir([BT_path 'CD_T_VP_' num2str(ij) '_' num2str(it) '.mat']);
        load([BT_path files.name],'p_res_change');
        for ip = 1:length(power_vec)
            for ipc = 1:length(p_change_vec)
                if capacity_idx == 1
                    p_response = zeros(length(Nvec_exp),length(mag_s));
                    for i = 1:length(Nvec_exp)                  
                        for jj = 1:length(mag_s)
%                             p_response(i,jj) = p_res_change(ip,ipc,max(Nvec_exp),i*2,jj);
                            p_response(i,jj) = p_res_change(ip,ipc,max(Nvec_exp),Nvec_exp(i),jj);
                        end    
                    end
                    p_no_response = 1-p_response;
                    p_response(p_response==0) = 1/(sampleSize_n*10);
                    p_no_response(p_no_response==0) = 1/(sampleSize_n*10);
                    p_log_like_response_yes = response_yes_trials.*log(p_response);
                    p_log_like_response_no = response_no_trials.*log(p_no_response);
                    p_log_like_mat = p_log_like_response_yes + p_log_like_response_no;
                    p_L(ij,it,ip,ipc) = sum(p_log_like_mat(:));
                elseif capacity_idx == 2    % fixed K
                    for ic = 1:length(capacity)
                        p_response = zeros(length(Nvec_exp),length(mag_s));
                        for i = 1:length(Nvec_exp)                  
                            for jj = 1:length(mag_s)
%                                 p_response(i,jj) = p_res_change(ip,ipc,ic,i*2,jj);
                                p_response(i,jj) = p_res_change(ip,ipc,ic,Nvec_exp(i),jj);
                            end    
                        end
                        p_no_response = 1-p_response;
                        p_response(p_response==0) = 1/(sampleSize_n*10);
                        p_no_response(p_no_response==0) = 1/(sampleSize_n*10);
                        p_log_like_response_yes = response_yes_trials.*log(p_response);
                        p_log_like_response_no = response_no_trials.*log(p_no_response);
                        p_log_like_mat = p_log_like_response_yes + p_log_like_response_no;
                        p_L(ij,it,ip,ipc,ic) = sum(p_log_like_mat(:));
                    end
                elseif capacity_idx == 3    % uniform K
                    for ic = 1:length(capacity)
                        p_response = zeros(length(Nvec_exp),length(mag_s));
                        Kmax = capacity(ic);
                        for i = 1:length(Nvec_exp)   
                            for jj = 1:length(mag_s)                                
                                for kk=1:Kmax
%                                     p_response(i,jj) = p_response(i,jj) + 1/Kmax*p_res_change(ip,ipc,kk,i*2,jj);
                                    p_response(i,jj) = p_response(i,jj) + 1/Kmax*p_res_change(ip,ipc,kk,Nvec_exp(i),jj);
                                end
                            end    
                        end
                        p_response = min(p_response,1);
                        p_no_response = 1-p_response;
                        p_response(p_response==0) = 1/(sampleSize_n*10);
                        p_no_response(p_no_response==0) = 1/(sampleSize_n*10);
                        p_log_like_response_yes = response_yes_trials.*log(p_response);
                        p_log_like_response_no = response_no_trials.*log(p_no_response);
                        p_log_like_mat = p_log_like_response_yes + p_log_like_response_no;
                        p_L(ij,it,ip,ipc,ic) = sum(p_log_like_mat(:));
                    end
                elseif capacity_idx == 4    % poisson K
                    for ic = 1:length(capacity_PK)
                        p_response = zeros(length(Nvec_exp),length(mag_s));
                        for i = 1:length(Nvec_exp)   
                            p_response_for_zero_K = 1/2;
%                             p_response_for_zero_K = 1/Nvec_exp(i);
                            for jj = 1:length(mag_s)
                                % case K=0
                                    p_response(i,jj) = p_response(i,jj) + poiss_w(ic,1)*p_response_for_zero_K;
                                % case K>0
                                for kk=1:8
%                                     p_response(i,jj) = p_response(i,jj) + poiss_w(ic,kk+1)*p_res_change(ip,ipc,kk,i*2,jj);
                                    p_response(i,jj) = p_response(i,jj) + poiss_w(ic,kk+1)*p_res_change(ip,ipc,kk,Nvec_exp(i),jj);
                                end
                            end    
                        end
                        p_response = min(p_response,1);
                        p_no_response = 1-p_response;
                        p_response(p_response==0) = 1/(sampleSize_n*10);
                        p_no_response(p_no_response==0) = 1/(sampleSize_n*10);
                        p_log_like_response_yes = response_yes_trials.*log(p_response);
                        p_log_like_response_no = response_no_trials.*log(p_no_response);
                        p_log_like_mat = p_log_like_response_yes + p_log_like_response_no;
                        p_L(ij,it,ip,ipc,ic) = sum(p_log_like_mat(:));
                    end
                end
            end
        end
    end
end

% model likelihood
p_L_max = max(p_L(:));
expL = exp(p_L-p_L_max);            % p landscape
if capacity_idx == 1
%     L = log(sum(expL(:))) + p_L_max - log(length(Jbar_vec)*length(tau_vec)*length(power_vec)*length(p_change_vec));
    L = log(trapz(Jbar_vec,trapz(tau_vec,trapz(power_vec,trapz(p_change_vec,expL,4),3),2),1)) + p_L_max - log(range(Jbar_vec)*range(tau_vec)*range(power_vec)*range(p_change_vec));
else
%     L = log(sum(expL(:))) + p_L_max - log(length(Jbar_vec)*length(tau_vec)*length(power_vec)*length(p_change_vec)*length(capacity_vec));
    L = log(trapz(Jbar_vec,trapz(tau_vec,trapz(power_vec,trapz(p_change_vec,trapz(capacity_vec,expL,5),4),3),2),1)) + p_L_max - log(range(Jbar_vec)*range(tau_vec)*range(power_vec)*range(p_change_vec)*range(capacity_vec));
end

% ML estimate
p_L = squeeze(p_L);
[M, I] = max(p_L(:));     

if capacity_idx == 1
    [idx_j, idx_t, idx_p, idx_pc] = ind2sub(size(p_L),I);
    ML_estimate = [Jbar_vec(idx_j) tau_vec(idx_t) power_vec(idx_p) p_change_vec(idx_pc) 8];
elseif capacity_idx > 1
    [idx_j, idx_t, idx_p, idx_pc, idx_c] = ind2sub(size(p_L),I);
    ML_estimate = [Jbar_vec(idx_j) tau_vec(idx_t) power_vec(idx_p) p_change_vec(idx_pc) capacity_vec(idx_c)];
end

% BIC
BIC = max(p_L(:)) - 3/2*log(length(data));
