function run_param_rec_test_de(nSbj)

PathName = 'de_delay/';

%% create fake data first

if ~exist('de_delay/de_delay_vp_test_1.mat','file')
    JbarRange = [0 100];
    PowerRange = [0 2];
    tauRange = [0 100];

    JbarVec = -min(JbarRange) + rand(1,nSbj)*range(JbarRange);
    PowerVec = -min(PowerRange) + rand(1,nSbj)*range(PowerRange);
    tauVec = -min(tauRange) + rand(1,nSbj)*range(tauRange);
    ntrials = 120;
%     ntrials = 1200;

    % interpolation of kappa over J
    upper_bound = 3000;
    k_vec = linspace(0,upper_bound,1e4);
    J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

    Nvec = [1 2 4 6];
    
    for i = 1:nSbj
        SbjIdx = i;
        J_bar = JbarVec(i);
        power = PowerVec(i);
        tau = tauVec(i);
        ParamInVec = [J_bar power tau];
        
        FileName = [PathName 'de_delay_vp_test_' num2str(SbjIdx) '.mat'];

        % pesudo random 
        SbjData = cell(1,4);
        for j = 1:length(Nvec)
            
            N = Nvec(j);
            J = gamrnd(J_bar*N^(-power)/tau, tau, 1, ntrials);
            kappa = interp1(J_vec,k_vec,J);
            kappa = min(kappa,700);

            temp = circ_vmrnd(0,kappa); % [-pi, pi]
            
            % remap
%             temp(temp > pi/2 & temp <= pi) = temp(temp > pi/2 & temp <= pi) - pi;
%             temp(temp > -pi & temp <= -pi/2) = temp(temp > -pi & temp <= -pi/2) + pi;
            
            SbjData{j} = temp;
            
        end

        save(FileName,'ParamInVec','SbjData');
    end
    
    fprintf('Fake subject files are created. \n');
    
%     for i = 1:nSbj
%         SbjIdx = i;
%         fit_de_delay_vp(SbjIdx)
%     end
else
    fprintf('Skipping fake subject generation... \n');
    
    ParamInVecMat = zeros(nSbj,3);
    ParamOutVecMat = zeros(nSbj,3);
    for i = 1:nSbj
        SbjIdx = i;
        FileName = [PathName 'T' num2str(SbjIdx) '-delay-errorfit' '.mat'];
        load(FileName,'ParamInVec','ParamOutVec');
        ParamInVecMat(i,:)= ParamInVec;
        ParamOutVecMat(i,:)= ParamOutVec;
    end
    
    figure;
    set(gca,'FontSize',11,'FontName','Arial');
    xfac = .13;
    yfac = .13;
    set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+1.6 yfac*5]);
    set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+1.6 yfac*5]);
    
    for i = 1:3
        subplot(1,3,i)
        plot(ParamInVecMat(:,i)',ParamOutVecMat(:,i)','o'); hold on
        plot(0:.1:max(ParamInVecMat(:,i)),0:.1:max(ParamInVecMat(:,i)),':');
    end
    
end


