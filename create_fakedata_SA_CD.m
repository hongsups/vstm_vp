function datamatrix = create_fakedata_SA_CD(ntrials,capacity_idx,J1,p_change,capacity,expid)

if capacity_idx == 1
    error('Capacity cannot be inifinite.');
end

% N and delta
% maxmag = pi;
maxmag = 90;
if expid == 1
    setsize_vec = randi(4,1,ntrials)*2;
else
    setsize_vec = 2.^(randi(4,1,ntrials)-1);
end    

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

% precomute K vectors
if capacity_idx == 2        % fixed K
    Kvec = repmat(capacity,1,ntrials);
elseif capacity_idx == 3    % uniform K
    Kvec = randi(capacity,1,ntrials);
elseif capacity_idx == 4    % poisson K
    Kvec = poissrnd(capacity,1,ntrials);
end

% determine changeness of trials in advance
change_idx_vec = randi(2,1,ntrials)-1;

datamatrix = zeros(ntrials, 6);
trialcnt = 0;
while (trialcnt<ntrials)
    trialcnt = trialcnt + 1;
    N = setsize_vec(trialcnt);                     % draw a random set size from Nvec    
   
    % choose target index and change magnitude
    change_idx = change_idx_vec(trialcnt);
    
    % choose target index and change magnitude
    target_idx = randi(N);      % target location
    mag = -maxmag + 2*maxmag*rand;
%     mag = maxmag*rand;          % change magnitude (uniform)
    mag = mag*change_idx;       % change present or not
    K = Kvec(trialcnt);
    
    if N >= K
        J1_mat = zeros(1, N);                   % initial: all zeros
        item_idx = 1:K;
        J1_mat(item_idx) = 1;                   % encode first K items with 1 chunk; rest doesnt receive chunk
    else
        J1_mat = repmat(floor(K/N),1,N);        % initial: all same resources
        item_idx = 1:mod(K, N);                 % pick first items where 1 will be added
        J1_mat(item_idx) = J1_mat(item_idx) + 1;% for those items, plus 1
    end    
   
    J = J1_mat*J1;

    k = interp1(J_vec,k_vec,J);
    dL_noise_mat = circ_vmrnd(0,k,2);
    dL_noise = (dL_noise_mat(:,1)-dL_noise_mat(:,2))';
    dL_noise(target_idx) = dL_noise(target_idx) + mag*pi/90;            % adding delta
    k_denom = sqrt(2*k.^2 + 2*k.^2.*cos(dL_noise));                     % trick to avoid NaN in bessel
    
    d_CL = besseli(0,k,1).^2.*exp(2*k-k_denom)./besseli(0,k_denom,1);
    d_CD = mean(d_CL)*p_change/(1-p_change);
    response_idx = d_CD>1;
    
    datamatrix(trialcnt,:) = [change_idx response_idx 0 0 N mag];  
end
