function CD_table_VP(Jbar_vec_idx,tau_vec_idx)

% parameter information
BT_path = 'CD/bigtable/';
fname = [BT_path 'CD_T_VP_' num2str(Jbar_vec_idx) '_' num2str(tau_vec_idx) '.mat'];

% if ~exist(fname,'file')

    Nvec = 1:8;                             % set size
    magstep = 1;                           % change magnitude step
    maxmag = 90;
    nsteps = 15;                            % # of steps in a parameter space
    sampleSize_n = 500;                       % sample size
    
%     magstep = 1;                           % change magnitude step
%     maxmag = 90;
%     nsteps = 3;                            % # of steps in a parameter space
%     sampleSize_n = 10;                       % sample size

    mag = 0:magstep:maxmag;                 % change magnitude
    % Jbar_vec = linspace(1, 80, nsteps);    
    % tau_vec = linspace(1, 60, nsteps);      
    Jbar_vec = logspace(log10(5),log10(300),nsteps);
    tau_vec = logspace(log10(5),log10(300),nsteps);
    power_vec = linspace(0, 2, nsteps);         
    p_change_vec = linspace(.2, .8, nsteps);    % prior
    capacity_vec = 1:max(Nvec);

    % interpolation of kappa over J
    upper_bound = 3000;
    k_vec = linspace(0,upper_bound,1e4);
    J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

    % ETL & preallocation
    totalIt = length(power_vec)*length(p_change_vec)*length(capacity_vec)*length(Nvec)*length(mag);
    p_res_change = zeros(length(power_vec),length(p_change_vec),length(capacity_vec),length(Nvec),length(mag));
    itCnt = 0;
    tic;

    for ij = Jbar_vec_idx %1:length(Jbar_vec)
        J1 = Jbar_vec(ij);
        for it = tau_vec_idx %1:length(tau_vec)
            t = tau_vec(it);
            for ip = 1:length(power_vec)
                power = power_vec(ip);
                for i = 1:length(Nvec)
                    N = Nvec(i);
                    for ic = 1:length(capacity_vec);
                        K = capacity_vec(ic);                                
                        J_bar = J1*min(K,N)^(-power);
                        for j = 1:length(mag)
                            delta = mag(j);
                            for ipc = 1:length(p_change_vec)
                                target_idx = randi(N,1,sampleSize_n);                                                                
                                p_change = p_change_vec(ipc);
                                J = gamrnd(J_bar/t, t, sampleSize_n, 2*N);              
                                k = interp1(J_vec,k_vec,J);
                                k1 = k(:, 1:N);
                                k2 = k(:, N+1:end);
                                
                                % x - y
                                dL_noise = circ_vmrnd(0,k1)-circ_vmrnd(0,k2);
                                for s=1:sampleSize_n
                                    dL_noise(s, target_idx(s)) = dL_noise(s, target_idx(s)) + delta*pi/90;
                                end
                                
                                % decision variable
                                k_denom = sqrt(k1.^2 + k2.^2 + 2*k1.*k2.*cos(dL_noise));
                                d_CL = besseli(0,k1,1).*besseli(0,k2,1).*exp(k1+k2-k_denom)./besseli(0,k_denom,1);
                                d_CL(:,1:max((N-K),0)) = 1;
                                d_CD = mean(d_CL',1)*p_change/(1-p_change);
                                decision_change = d_CD>1;
                                
                                p_res_change(ip, ipc, ic, i, j) = mean(decision_change);    
                                
                                itCnt=itCnt+1;
                            end
                        end
                        fprintf(['table_VP' ' ETL=%2.1f min\n'],toc/itCnt*(totalIt-itCnt)/60);
                    end
                end
            end
        end
    end

    % save
    save(fname);
% else
%     fprintf('File already exists - skipping\n');
% end

    