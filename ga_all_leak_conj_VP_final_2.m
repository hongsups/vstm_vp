function ga_all_leak_conj_VP_final_2(test_idx,sbj_idx,omega_idx,power_idx,popSize,genNum,sampleSize,ntrials)

if ~exist('omega_idx','var')
    omega_idx = 0;
end
if ~exist('power_idx','var')
    power_idx = 0;
end

model_idx = omega_idx*10 + power_idx;
% o = 0, p = 0 -> 0
% o = 1, p = 0 -> 10
% o = 2, p = 0 -> 20
% o = 0, p = 1 -> 1
% o = 1, p = 1 -> 11
% o = 2, p = 1 -> 21
% o = 0, p = 2 -> 2
% o = 1, p = 2 -> 12
% o = 2, p = 2 -> 22
% omega 1, omega 2, power 1, power 2

% Here you can both compute LLH for shared.
% The idea is, you first draw Js from Jbar of feature 1 and 2.
% J with nonzero delta: this is for one-feature condition.
% For two-feature condition, we need a zero delta for the irrelevant.

% 5 parameters for shared: Jbar1, tau1, Jbar2, tau2, rho (all vectors are in this order)
% 4 parameters for independent

% load subject data here and use it as a passing input for ga
if test_idx == 1 % fake data test
    fname = ['ga_conj_VP_test_all_leak_o' num2str(omega_idx) '_p' num2str(power_idx) '_' num2str(sbj_idx) '_' num2str(popSize) '_' num2str(genNum) '_' num2str(sampleSize) '_' num2str(ntrials) '.mat'];
    sbj_name = ['test_conj_VP_o' num2str(omega_idx) '_p' num2str(power_idx) '_' num2str(sbj_idx) '_' num2str(popSize) '_' num2str(genNum) '_' num2str(sampleSize) '_' num2str(ntrials)];
    load(fname);
    data4_sbj = data4;
    data8_sbj = data8;
    dataS_sbj = dataS;
elseif test_idx == 0
   % load subject data
   subjects = {'AAR', 'HS0', 'WK0', 'HJK', 'JSK'};
   subjid = subjects{sbj_idx};
   sbj_name = ['ga_conj_VP_o' num2str(omega_idx) '_p' num2str(power_idx) '_' subjid '_' num2str(popSize) '_' num2str(genNum) '_' num2str(sampleSize) '_' num2str(ntrials)];
   [data4, data8, dataS] = loaddata_conjAll(subjid);
    data4_sbj = data4;
    data8_sbj = data8;
    dataS_sbj = dataS;
end

% set population size (default: 20, but with small # of params, 10 is enough)
% set initial range (default is [0 1]);
% Handle to the function that produces mutation children
% Mutation function: default (Gaussian)
% StallGenLimit = ## (after ## gen of stall, it stops)

% popSize = 64;
% genNum = 64;
n_add_param = omega_idx + power_idx;

LB = [1 1 1 1 zeros(1,n_add_param)];
UB = [200 200 200 200 ones(1,omega_idx) ones(1,power_idx)*2];
initRange = [LB; UB];
nvar = 4 + n_add_param;
    
% plot the process
% opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping});
% set popSize, initial range, stallgenlim, gennum
opts = gaoptimset('PopulationSize',popSize,'PopInitRange',initRange,'Generations',genNum);

% run ga solver: ub, lb, passing subject index as an additional variable
FitnessFun = {@LLH,data4_sbj,data8_sbj,dataS_sbj,model_idx,sampleSize};
[x, fval, ~, Output] = ga(FitnessFun,nvar,[],[],[],[],LB,UB,[],[],opts);

fprintf('The number of generations was : %d\n', Output.generations);
fprintf('The number of function evaluations was : %d\n', Output.funccount);
fprintf('The best function value found was : %g\n', -fval);

% LLH input: Jbar1, tau1, Jbar2, tau2, omeag1, omeag2, power1, power2)

% MLE (output)
Jbar1_out = x(1);
tau1_out = x(2);
Jbar2_out = x(3);
tau2_out = x(4);

if omega_idx == 0 && power_idx == 0
%     par_rest_vec 
    omega1_out = 0;
    omega2_out = 0;
    
    power1_out = 1;
    power2_out = 1;
else
    % based on param_vec, assign proper values to parameters
    
    par_rest_vec = x(5:end);
    omega_vec = par_rest_vec(1:omega_idx);
    power_vec = par_rest_vec(omega_idx+1:end);

    if omega_idx == 0
        omega1_out = 0;
        omega2_out = 0;
    elseif length(omega_vec) == 1
        omega1_out = omega_vec;
        omega2_out = omega_vec;
    else
        omega1_out = omega_vec(1);
        omega2_out = omega_vec(2);
    end

    if power_idx == 0
        power1_out = 1;
        power2_out = 1;
    elseif length(power_vec) == 1
        power1_out = power_vec;
        power2_out = power_vec;
    else
        power1_out = power_vec(1);
        power2_out = power_vec(2);
    end
    
end


param_out_vec = x(1:nvar);
if test_idx == 1
    param_in_vec = [Jbar1_in tau1_in Jbar2_in tau2_in omega1_in omega2_in power1_in power2_in];
elseif test_idx == 0
    param_in_vec = [];
end


% generate fake data
paramvec = [Jbar1_out,Jbar2_out,tau1_out,tau2_out,power1_out,power2_out];

fname_temp = ['temp_' sbj_name  '.mat'];
save(fname_temp);

data4_model = create_fakedata_V1_IR_conj(ntrials,paramvec,4);
data8_model = create_fakedata_V1_IR_conj(ntrials,paramvec,8);
dataS_model = create_fakedata_VPsplitleak(ntrials,[Jbar1_out,Jbar2_out,tau1_out,tau2_out,omega1_out,omega2_out,power1_out,power2_out]);

generate_modelfit_VP_all_leak(sbj_name,data4_sbj,data8_sbj,dataS_sbj,data4_model,data8_model,dataS_model,param_in_vec,param_out_vec,fval);

function LLH_ga_VP = LLH(param_vec,data4,data8,dataS,model_idx,sampleSize)

Jbar1_out = param_vec(1);
tau1_out = param_vec(2);
Jbar2_out = param_vec(3);
tau2_out = param_vec(4);
if length(param_vec) > 4
    param_rest_vec = param_vec(5:end);
end
% 1st digit of the model_idx = # of omega (modular of 10)
% 2nd digit of the model_idx = # of omega 
% param order in param_vec: Jbar1, tau1, Jbar2, tau2, omega1, omega2, power1, power2

% first count the number of additional parameters (omegas, powers)
num_omega = floor(model_idx/10);
num_power = model_idx - num_omega*10;

% error warning
if length(param_vec)-4 ~= num_omega+num_power
    save alles
    error('Parameter numbers do not match.');
end

% based on param_vec, assign proper values to parameters
if num_omega == 0 && num_power == 0
    
else
    omega_vec = param_rest_vec(1:num_omega);
    power_vec = param_rest_vec(num_omega+1:end);
end

if num_omega == 0
    omega1_out = 0;
    omega2_out = 0;
elseif length(omega_vec) == 1
    omega1_out = omega_vec;
    omega2_out = omega_vec;
else
    omega1_out = omega_vec(1);
    omega2_out = omega_vec(2);
end

if num_power == 0
    power1_out = 1;
    power2_out = 1;
elseif length(power_vec) == 1
    power1_out = power_vec;
    power2_out = power_vec;
else
    power1_out = power_vec(1);
    power2_out = power_vec(2);
end

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

% sample size for MC simulation
delta_vec = linspace(0, 90, 91);

% preallocation
target_idx_conj = 1;
p_corr_S = zeros(1,length(delta_vec));
dmat1_4 = zeros(length(delta_vec),sampleSize,4);
dmat2_4 = zeros(length(delta_vec),sampleSize,4);
dmat1_8 = zeros(length(delta_vec),sampleSize,8);
dmat2_8 = zeros(length(delta_vec),sampleSize,8);
for m = 1:length(delta_vec)
    delta = delta_vec(m);

    %% N = 4, two feature
    % color: odd number, orientation: even number
    N = 4;
    J1 = gamrnd(Jbar1_out*N^(-power1_out)/tau1_out, tau1_out, sampleSize, 2*N);
    J2 = gamrnd(Jbar2_out*N^(-power2_out)/tau2_out, tau2_out, sampleSize, 2*N);
    k1 = interp1(J_vec,k_vec,J1);
    k11 = k1(:, 1:N);
    k12 = k1(:, N+1:end);
    k2 = interp1(J_vec,k_vec,J2);
    k21 = k2(:, 1:N);
    k22 = k2(:, N+1:end);
    dL_noise1 = circ_vmrnd(0,k11)-circ_vmrnd(0,k12);    % size: nsample x N
    dL_noise2 = circ_vmrnd(0,k21)-circ_vmrnd(0,k22);
    
    dL_noise1(:, target_idx_conj) = dL_noise1(:, target_idx_conj) + delta*pi/90;
    dL_noise2(:, target_idx_conj) = dL_noise2(:, target_idx_conj) + delta*pi/90;
    k_denom1 = sqrt(k11.^2 + k12.^2 + 2*k11.*k12.*cos(dL_noise1));
    k_denom2 = sqrt(k21.^2 + k22.^2 + 2*k21.*k22.*cos(dL_noise2));
    d1 = besseli(0,k11,1).*besseli(0,k12,1).*exp(k11+k12-k_denom1)./besseli(0,k_denom1,1);
    d2 = besseli(0,k21,1).*besseli(0,k22,1).*exp(k21+k22-k_denom2)./besseli(0,k_denom2,1);
    dmat1_4(m,:,:) = d1;
    dmat2_4(m,:,:) = d2;

    %% Split and leak
    J1 = gamrnd(Jbar1_out*N^(-power1_out)*(1-omega1_out)/tau1_out, tau1_out, sampleSize, 2*N);
    J2 = gamrnd(Jbar2_out*N^(-power2_out)*(1-omega2_out)/tau2_out, tau2_out, sampleSize, 2*N);
    k1 = interp1(J_vec,k_vec,J1);
    k11 = k1(:, 1:N);
    k12 = k1(:, N+1:end);
    k2 = interp1(J_vec,k_vec,J2);
    k21 = k2(:, 1:N);
    k22 = k2(:, N+1:end);
    dL_noise1 = circ_vmrnd(0,k11)-circ_vmrnd(0,k12);    % size: nsample x N
    dL_noise2 = circ_vmrnd(0,k21)-circ_vmrnd(0,k22);

    dL_noise = [dL_noise1(:,1) dL_noise2(:,1) dL_noise1(:,2) dL_noise2(:,2) dL_noise1(:,3) dL_noise2(:,3) dL_noise1(:,4) dL_noise2(:,4)];                   % combine: nsample x 2N (8)    
    target_idx = randi(2*N,1,sampleSize);
    for s=1:sampleSize
        dL_noise(s, target_idx(s)) = dL_noise(s, target_idx(s)) + delta*pi/90;
    end
    dL_noise_1 = dL_noise(:, [1 3 5 7]);
    dL_noise_2 = dL_noise(:, [2 4 6 8]);
    k_denom1 = sqrt(k11.^2 + k12.^2 + 2*k11.*k12.*cos(dL_noise_1));                            
    k_denom2 = sqrt(k21.^2 + k22.^2 + 2*k21.*k22.*cos(dL_noise_2));                            
    d1 = besseli(0,k11,1).*besseli(0,k12,1).*exp(k11+k12-k_denom1)./besseli(0,k_denom1,1);
    d2 = besseli(0,k21,1).*besseli(0,k22,1).*exp(k21+k22-k_denom2)./besseli(0,k_denom2,1);
    d = [d1(:,1) d2(:,1) d1(:,2) d2(:,2) d1(:,3) d2(:,3) d1(:,4) d2(:,4)];    
    c = zeros(1,sampleSize);
    for s=1:sampleSize
        max_idx = find(d(s,:)==max(d(s,:)));
        c(s) = target_idx(s) == max_idx(randi(length(max_idx)));  %d(target_idx)==max(d);
    end
    p_corr_S(m) = mean(c);    
    
    %% N = 8, two feature
    N = 8;
    J1 = gamrnd(Jbar1_out*N^(-power1_out)/tau1_out, tau1_out, sampleSize, 2*N);
    J2 = gamrnd(Jbar2_out*N^(-power2_out)/tau2_out, tau2_out, sampleSize, 2*N);
    k1 = interp1(J_vec,k_vec,J1);
    k11 = k1(:, 1:N);
    k12 = k1(:, N+1:end);
    k2 = interp1(J_vec,k_vec,J2);
    k21 = k2(:, 1:N);
    k22 = k2(:, N+1:end);
    dL_noise1 = circ_vmrnd(0,k11)-circ_vmrnd(0,k12);    % size: nsample x N
    dL_noise2 = circ_vmrnd(0,k21)-circ_vmrnd(0,k22);
    dL_noise1(:, target_idx_conj) = dL_noise1(:, target_idx_conj) + delta*pi/90;
    dL_noise2(:, target_idx_conj) = dL_noise2(:, target_idx_conj) + delta*pi/90;
    k_denom1 = sqrt(k11.^2 + k12.^2 + 2*k11.*k12.*cos(dL_noise1));
    k_denom2 = sqrt(k21.^2 + k22.^2 + 2*k21.*k22.*cos(dL_noise2));
    d1 = besseli(0,k11,1).*besseli(0,k12,1).*exp(k11+k12-k_denom1)./besseli(0,k_denom1,1);
    d2 = besseli(0,k21,1).*besseli(0,k22,1).*exp(k21+k22-k_denom2)./besseli(0,k_denom2,1);
    dmat1_8(m,:,:) = d1;
    dmat2_8(m,:,:) = d2;
    
end
% TWO-FEATURE predictions for FEATURE 1 (Feature 2 = zero)
% replicate zero-delta for all non-zero deltas
dmat2_zero = repmat(dmat2_4(1,:,:),length(delta_vec),1);  
% compute the new decision variable (mean)
d_conj4 = (dmat1_4 + dmat2_zero)/2;
% add noise for convenience
d_conj4 = normrnd(d_conj4,1e-5);
% apply decision rule
c_conj4 = d_conj4(:,:,target_idx_conj)==max(d_conj4,[],3);
% compute mean to calculate the predicion (over samples)
p_corr_two1_4 = mean(c_conj4');

% TWO-FEATURE predictions for FEATURE 2 (Feature 1 = zero)
% replicate zero-delta for all non-zero deltas
dmat1_zero = repmat(dmat1_4(1,:,:),length(delta_vec),1);  
% compute the new decision variable (mean)
d_conj4 = (dmat2_4 + dmat1_zero)/2;
% add noise for convenience
d_conj4 = normrnd(d_conj4,1e-5);
% apply decision rule
c_conj4 = d_conj4(:,:,target_idx_conj)==max(d_conj4,[],3);
% compute mean to calculate the predicion (over samples)
p_corr_two2_4 = mean(c_conj4');

% TWO-FEATURE predictions for FEATURE 1 (Feature 2 = zero)
% replicate zero-delta for all non-zero deltas
dmat2_zero = repmat(dmat2_8(1,:,:),length(delta_vec),1);  
% compute the new decision variable (mean)
d_conj8 = (dmat1_8 + dmat2_zero)/2;
% add noise for convenience
d_conj8 = normrnd(d_conj8,1e-5);
% apply decision rule
c_conj8 = d_conj8(:,:,target_idx_conj)==max(d_conj8,[],3);
% compute mean to calculate the predicion (over samples)
p_corr_two1_8 = mean(c_conj8');

% TWO-FEATURE predictions for FEATURE 2 (Feature 1 = zero)
% replicate zero-delta for all non-zero deltas
dmat1_zero = repmat(dmat1_8(1,:,:),length(delta_vec),1);  
% compute the new decision variable (mean)
d_conj8 = (dmat2_8 + dmat1_zero)/2;
% add noise for convenience
d_conj8 = normrnd(d_conj8,1e-5);
% apply decision rule
c_conj8 = d_conj8(:,:,target_idx_conj)==max(d_conj8,[],3);
% compute mean to calculate the predicion (over samples)
p_corr_two2_8 = mean(c_conj8');

% load subject data: 4
conj_idx_ori = data4(:,2)~=0;  % indices of conjunction trials in which change was in 1st dimension
conj_idx_col = data4(:,3)~=0;  % indices of conjunction trials in which change was in 2nd dimension
magidx_conj_ori = round(abs(data4(conj_idx_ori,2)));
magidx_conj_col = round(abs(data4(conj_idx_col,3)));
perfidx_conj_ori = data4(conj_idx_ori, 9);
perfidx_conj_col = data4(conj_idx_col, 9);
corr_trials_conj_ori = zeros(1,length(delta_vec));
incorr_trials_conj_ori = zeros(1,length(delta_vec));
corr_trials_conj_col = zeros(1,length(delta_vec));
incorr_trials_conj_col = zeros(1,length(delta_vec));
for i = 1:length(delta_vec)
    idx_conj_ori = find(magidx_conj_ori == delta_vec(i));
    idx_conj_col = find(magidx_conj_col == delta_vec(i));
    corr_trials_conj_ori(i) = sum(perfidx_conj_ori(idx_conj_ori));
    incorr_trials_conj_ori(i) = sum(perfidx_conj_ori(idx_conj_ori)==0);
    corr_trials_conj_col(i) = sum(perfidx_conj_col(idx_conj_col));
    incorr_trials_conj_col(i) = sum(perfidx_conj_col(idx_conj_col)==0);
end
p_incorr_two1_4 = 1-p_corr_two1_4;
p_corr_two1_4(p_corr_two1_4==0) = 1/(sampleSize*10);
p_incorr_two1_4(p_incorr_two1_4==0) = 1/(sampleSize*10);
p_log_like_corr_two1_4 = corr_trials_conj_ori.*log(p_corr_two1_4);
p_log_like_incorr_two1_4 = incorr_trials_conj_ori.*log(p_incorr_two1_4);
p_log_like_conj_ori = p_log_like_corr_two1_4 + p_log_like_incorr_two1_4;
p_incorr_two2_4 = 1-p_corr_two2_4;
p_corr_two2_4(p_corr_two2_4==0) = 1/(sampleSize*10);
p_incorr_two2_4(p_incorr_two2_4==0) = 1/(sampleSize*10);
p_log_like_corr_two2_4 = corr_trials_conj_col.*log(p_corr_two2_4);
p_log_like_incorr_two2_4 = incorr_trials_conj_col.*log(p_incorr_two2_4);
p_log_like_conj_col = p_log_like_corr_two2_4 + p_log_like_incorr_two2_4;
p_LLH_4 = sum(p_log_like_conj_ori(:)) + sum(p_log_like_conj_col(:));

% load subject data: 8
conj_idx_ori = data8(:,2)~=0;  % indices of conjunction trials in which change was in 1st dimension
conj_idx_col = data8(:,3)~=0;  % indices of conjunction trials in which change was in 2nd dimension
magidx_conj_ori = round(abs(data8(conj_idx_ori,2)));
magidx_conj_col = round(abs(data8(conj_idx_col,3)));
perfidx_conj_ori = data8(conj_idx_ori, 9);
perfidx_conj_col = data8(conj_idx_col, 9);
corr_trials_conj_ori = zeros(1,length(delta_vec));
incorr_trials_conj_ori = zeros(1,length(delta_vec));
corr_trials_conj_col = zeros(1,length(delta_vec));
incorr_trials_conj_col = zeros(1,length(delta_vec));
for i = 1:length(delta_vec)
    idx_conj_ori = find(magidx_conj_ori == delta_vec(i));
    idx_conj_col = find(magidx_conj_col == delta_vec(i));
    corr_trials_conj_ori(i) = sum(perfidx_conj_ori(idx_conj_ori));
    incorr_trials_conj_ori(i) = sum(perfidx_conj_ori(idx_conj_ori)==0);
    corr_trials_conj_col(i) = sum(perfidx_conj_col(idx_conj_col));
    incorr_trials_conj_col(i) = sum(perfidx_conj_col(idx_conj_col)==0);
end
p_incorr_two1_8 = 1-p_corr_two1_8;
p_corr_two1_8(p_corr_two1_8==0) = 1/(sampleSize*10);
p_incorr_two1_8(p_incorr_two1_8==0) = 1/(sampleSize*10);
p_log_like_corr_two1_8 = corr_trials_conj_ori.*log(p_corr_two1_8);
p_log_like_incorr_two1_8 = incorr_trials_conj_ori.*log(p_incorr_two1_8);
p_log_like_conj_ori = p_log_like_corr_two1_8 + p_log_like_incorr_two1_8;
p_incorr_two2_8 = 1-p_corr_two2_8;
p_corr_two2_8(p_corr_two2_8==0) = 1/(sampleSize*10);
p_incorr_two2_8(p_incorr_two2_8==0) = 1/(sampleSize*10);
p_log_like_corr_two2_8 = corr_trials_conj_col.*log(p_corr_two2_8);
p_log_like_incorr_two2_8 = incorr_trials_conj_col.*log(p_incorr_two2_8);
p_log_like_conj_col = p_log_like_corr_two2_8 + p_log_like_incorr_two2_8;
p_LLH_8 = sum(p_log_like_conj_ori(:)) + sum(p_log_like_conj_col(:));

% load subject data: S
data_delta = dataS(:,2) + dataS(:,3);
magidx = round(abs(data_delta));   
perfidx = dataS(:, 9);
corr_trials = zeros(1,length(delta_vec));
incorr_trials = zeros(1,length(delta_vec));
for jj = 1:length(delta_vec)
    idx = find(magidx == delta_vec(jj));
    corr_trials(jj) = sum(perfidx(idx));               % # of corr. trials
    incorr_trials(jj) = sum(perfidx(idx)==0);          % # of incorr. trials
end    
p_incorr_S = 1-p_corr_S;
p_corr_S(p_corr_S==0) = 1/(sampleSize*10);
p_incorr_S(p_incorr_S==0) = 1/(sampleSize*10);
p_log_like_corr_S = corr_trials.*log(p_corr_S);
p_log_like_incorr_S = incorr_trials.*log(p_incorr_S);
p_log_like_mat_S = p_log_like_corr_S + p_log_like_incorr_S;
p_LLH_S = sum(p_log_like_mat_S(:));

% since GA FINDS MINIMUM, we need to flipt the sign.
LLH_ga_VP = -(p_LLH_4 + p_LLH_8 + p_LLH_S);

function generate_modelfit_VP_all_leak(sbj_name,data4_sbj,data8_sbj,dataS_sbj,data4_model,data8_model,dataS_model,param_in_vec,param_out_vec,fval)

nBins = 9;
[bin, ~] = binning(nBins);

% N = 4
perf4_sbj = performance_conj(data4_sbj,8,nBins);             % has both orientatio and color!
perf4_model = performance_conj(data4_model,8,nBins);    % has both orientatio and color!

perf4_sbj_ori = perf4_sbj(1,:);
perf4_sbj_col = perf4_sbj(2,:);
perf4_model_ori = perf4_model(1,:);
perf4_model_col = perf4_model(2,:);

% N = 8
perf8_sbj = performance_conj(data8_sbj,8,nBins);             % has both orientatio and color!
perf8_model = performance_conj(data8_model,8,nBins);    % has both orientatio and color!

perf8_sbj_ori = perf8_sbj(1,:);
perf8_sbj_col = perf8_sbj(2,:);
perf8_model_ori = perf8_model(1,:);
perf8_model_col = perf8_model(2,:);

% Split
dataS_sbj_ori = dataS_sbj(dataS_sbj(:,3)==0,:);
dataS_sbj_col = dataS_sbj(dataS_sbj(:,2)==0,:);
dataS_model_ori = dataS_model(dataS_model(:,3)==0,:);
dataS_model_col = dataS_model(dataS_model(:,2)==0,:);

perfS_sbj_ori = performance_conj(dataS_sbj_ori,8,nBins);    
perfS_sbj_col = performance_conj(dataS_sbj_col,8,nBins);    
perfS_sbj_ori = perfS_sbj_ori(1,:);
perfS_sbj_col = perfS_sbj_col(2,:);

perfS_model_ori = performance_conj(dataS_model_ori,8,nBins);
perfS_model_col = performance_conj(dataS_model_col,8,nBins);    
perfS_model_ori = perfS_model_ori(1,:);
perfS_model_col = perfS_model_col(2,:);

figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1.2;
yfac = .2;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+.5 yfac*5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+.5 yfac*5]);

subplot(1,2,1);     % orientation
plot(bin,perf4_sbj_ori,'bo');hold on
plot(bin,perf8_sbj_ori,'ro');hold on
plot(bin,perfS_sbj_ori,'go');hold on
plot(bin,perf4_model_ori,'b');hold on
plot(bin,perf8_model_ori,'r');hold on
plot(bin,perfS_model_ori,'g');
title('Orientation');xlabel('Change magnitude ({\circ})');ylabel('Proportion correct');xlim([0 90]);ylim([0 1]);hold on;
set(gca,'YTick',0:.2:1)
set(gca,'XTick',0:30:90)

subplot(1,2,2);     % color
plot(bin,perf4_sbj_col,'bo');hold on
plot(bin,perf8_sbj_col,'ro');hold on
plot(bin,perfS_sbj_col,'go');hold on
plot(bin,perf4_model_col,'b');hold on
plot(bin,perf8_model_col,'r');hold on
plot(bin,perfS_model_col,'g');
title('Color');xlabel('Change magnitude ({\circ})');ylabel('Proportion correct');xlim([0 90]);ylim([0 1]);hold on;
set(gca,'YTick',0:.2:1)
set(gca,'XTick',0:30:90)

if isempty(param_in_vec)==0
    text(1, .95, ['IN: ' num2str(param_in_vec, 3)]);
end
text(1, .90, ['OUT: ' num2str(param_out_vec, 3)]);
saveas(gcf,['fit_' sbj_name '.fig'],'fig');
saveas(gcf,['fit_' sbj_name '.emf'],'emf');

fname = [sbj_name '.mat'];
save(fname);
