function [L, ML_estimate, BIC] = create_LLH_SA_CD(capacity_idx,expid,subjid)
    
if capacity_idx == 1
    error('Capacity cannot be inifinite.');
end

% load parameters from the table file
BT_path = 'CD/bigtable/';
files = dir([BT_path 'CD_T_SA_1.mat']);
load([BT_path files.name],'nsteps','J1_vec','p_change_vec','magstep','mag','sampleSize_n','Nvec','capacity_SA');

if expid == 1
    Nvec_exp = 2:2:8;                       % set sizes used in the experiments
else
    Nvec_exp = [1 2 4 8];
end
capacity_PK = linspace(1/3,10,nsteps);  % PK

if capacity_idx == 4
    capacity_vec = capacity_PK;
else
    capacity_vec = capacity_SA;
end

% load subject data
data = loaddata_CD(subjid,expid);
if expid < 3
    change_idx_vec = ceil((data(:,1))/max(data(:,1)));
    mag_vec = data(:,6)/pi*180.*change_idx_vec;     % [0,180] (You need to multiply change_idx)    
    mag_vec(mag_vec>90)=180-mag_vec(mag_vec>90);    % ceil(mag_vec) -> min: 0, max: 90
else
    mag_vec = abs(data(:,6));    % For fake data, you don't have to multiply change_idx_vec. It's already take care of in fake data generatio code.
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
poiss_w = zeros(length(capacity_PK),maxKinPoisson+1);
for ii=1:length(capacity_PK)
    poiss_w(ii,:) = poisspdf(0:maxKinPoisson,capacity_PK(ii));
    poiss_w(ii,:) = poiss_w(ii,:)/sum(poiss_w(ii,:));
end

% p_res_change(ipc, ic, i, j)
% compute likelihood
p_L = zeros(length(J1_vec),length(p_change_vec),length(capacity_vec));
for ij = 1:length(J1_vec)
    % load a lookup table
    files = dir([BT_path 'CD_T_SA_' num2str(ij)  '.mat']);
    load([BT_path files.name],'p_res_change');
    for ipc = 1:length(p_change_vec)
        if capacity_idx == 2                    % fixed K
            for ic = 1:length(capacity_SA)                
                p_response = zeros(length(Nvec_exp),length(mag_s));
                for i = 1:length(Nvec_exp)
                    for jj = 1:length(mag_s)                        
%                         p_response(i,jj) = p_res_change(ipc,ic,i*2,jj);
                        p_response(i,jj) = p_res_change(ipc,ic,Nvec_exp(i),jj);                        
                    end
                end
                
                p_no_response = 1-p_response;
                p_response(p_response==0) = 1/(sampleSize_n*10);
                p_no_response(p_no_response==0) = 1/(sampleSize_n*10);

                p_log_like_response_yes = response_yes_trials.*log(p_response);
                p_log_like_response_no = response_no_trials.*log(p_no_response);
                p_log_like_mat = p_log_like_response_yes + p_log_like_response_no;
                p_L(ij,ipc,ic) = sum(p_log_like_mat(:));
            end
        elseif capacity_idx == 3                % uniform K
            for ic = 1:length(capacity_SA)
                p_response = zeros(length(Nvec_exp),length(mag_s));
                Kmax = capacity_SA(ic);
                for i = 1:length(Nvec_exp)
                    for jj = 1:length(mag_s)
                        for kk=1:Kmax
%                             nEncoded = min(kk,Nvec_exp(i));  % number of encoded items
%                             p_response(i,jj) = p_response(i,jj) + 1/Kmax*p_res_change(ipc,nEncoded,i*2,jj);
%                             p_response(i,jj) = p_response(i,jj) + 1/Kmax*p_res_change(ipc,kk,i*2,jj);
                            p_response(i,jj) = p_response(i,jj) + 1/Kmax*p_res_change(ipc,kk,Nvec_exp(i),jj);                            
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
                p_L(ij,ipc,ic) = sum(p_log_like_mat(:));
            end
        elseif capacity_idx == 4                  % poisson K
            for ic = 1:length(capacity_PK)
                p_response = zeros(length(Nvec_exp),length(mag_s));
                for i = 1:length(Nvec_exp)
%                     p_response_for_zero_K = 1/Nvec_exp(i);
                    p_response_for_zero_K = 1/2;                    
                    for jj = 1:length(mag_s)
                        % case K=0
                        p_response(i,jj) = p_response(i,jj) + poiss_w(ic,1)*p_response_for_zero_K;
                        % case K>0
                        for kk=1:maxKinPoisson
%                             p_response(i,jj) = p_response(i,jj) + poiss_w(ic,kk+1)*p_res_change(ipc,kk,i*2,jj);
                            p_response(i,jj) = p_response(i,jj) + poiss_w(ic,kk+1)*p_res_change(ipc,kk,Nvec_exp(i),jj);
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
                p_L(ij,ipc,ic) = sum(p_log_like_mat(:));
            end
        end
    end
end

% model likelihood
p_L_max = max(p_L(:));
expL = exp(p_L-p_L_max);            % p landscape

% linear
% L = log(sum(expL(:))) + p_L_max - log(length(J1_vec)*length(capacity_vec)*length(p_change_vec));

% log
L = log(trapz(J1_vec,trapz(p_change_vec,trapz(capacity_vec,expL,3),2),1)) + p_L_max - log(range(J1_vec)*range(capacity_vec)*range(p_change_vec));

% ML estimate
p_L = squeeze(p_L);
[M, I] = max(p_L(:));     
[idx_j, idx_pc, idx_c] = ind2sub(size(p_L),I);
ML_estimate = [J1_vec(idx_j) p_change_vec(idx_pc) capacity_vec(idx_c)];

% BIC
BIC = max(p_L(:)) - 3/2*log(length(data)); 