function datamatrix = create_fakedata_VP_CD(ntrials,capacity_idx,J_bar,t,power,p_change,capacity,expid)

% set size and delta vectors
% maxmag = pi;
% N and delta
maxmag = 90;
if expid == 1
    setsize_vec = randi(4,1,ntrials)*2;
    Nvec_exp = 2:2:8;
else
    setsize_vec = 2.^(randi(4,1,ntrials)-1);
    Nvec_exp = [1 2 4 8];
end    

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

% precomute K vectors
if capacity_idx == 1        % infinite K
    Kvec = repmat(max(Nvec_exp),1,ntrials);
elseif capacity_idx == 2    % fixed K
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

    J = gamrnd((J_bar*min(K,N)^(-power))/t, t, 1, 2*N);                            
    k = interp1(J_vec,k_vec,J);
    k1 = k(1:N);
    k2 = k(N+1:end);
    dL_noise = circ_vmrnd(0,k1)-circ_vmrnd(0,k2);
    dL_noise(target_idx) = dL_noise(target_idx) + mag*pi/90;

    k_denom = sqrt(k1.^2 + k2.^2 + 2*k1.*k2.*cos(dL_noise));
    d_CL = besseli(0,k1,1).*besseli(0,k2,1).*exp(k1+k2-k_denom)./besseli(0,k_denom,1);
    d_CL(1:max((N-K),0)) = 1;
    
    d_CD = mean(d_CL)*p_change/(1-p_change);
    response_idx = d_CD>1;

    datamatrix(trialcnt,:) = [change_idx response_idx 0 0 N mag];  
end
