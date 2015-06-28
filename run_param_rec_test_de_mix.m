function run_param_rec_test_de_mix(nSbj)

% need to be implemented

PathName = 'de_delay/';

%% create fake data first

if ~exist('de_delay/de_delay_mix_test_1.mat','file')
    JRange = [0 100];
    wRange = [0 1];

    JVec = -min(JRange) + rand(1,nSbj)*range(JRange);
    wVec = -min(wRange) + rand(1,nSbj)*range(wRange);
    ntrials = 120;

    % interpolation of kappa over J
    upper_bound = 3000;
    k_vec = linspace(0,upper_bound,1e4);
    J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);
    Nvec = [1 2 4 6];

    power = 1;
    for i = 1:nSbj
        SbjIdx = i;
        J1 = JVec(i);
        w = wVec(i);
        
        ParamInVec = [J1 w];
        FileName = [PathName 'de_delay_mix_test_' num2str(SbjIdx) '.mat'];
        
        SbjData = cell(1,4);
        for j = 1:length(Nvec)
            N = Nvec(j);
            J = J1*N^(-power);
            kappa = interp1(J_vec,k_vec,J);
            kappa = min(kappa,700);
                
            SbjDataVec = zeros(1,ntrials);
            for k = 1:ntrials
                randidx = rand;
                if randidx <= w % VM
                    SbjDataVec(k) = circ_vmrnd(0,kappa,1);
                else
                    uniformrange = [-pi pi];
                    SbjDataVec(k) = min(uniformrange) + rand*range(uniformrange);
                end
            end
            SbjData{j} = SbjDataVec;
        end

        save(FileName,'ParamInVec','SbjData');
    end
    
    fprintf('Fake subject files are created. \n');

%     for i = 1:nSbj
%         SbjIdx = i;
%         fit_de_delay_mix(SbjIdx)
%     end
    
else
    fprintf('Skipping fake subject generation... \n');
    
    ParamInVecMat = zeros(nSbj,2);
    ParamOutVecMat = zeros(nSbj,2);
    for i = 1:nSbj
        SbjIdx = i;
        FileName = [PathName 'T' num2str(SbjIdx) '-delay-mix-errorfit' '.mat'];
        load(FileName,'ParamInVec','ParamOutVec');
        
        ParamInVecMat(i,:)= ParamInVec;
        ParamOutVecMat(i,:)= ParamOutVec;
    end
    
    figure;
    set(gca,'FontSize',11,'FontName','Arial');
    xfac = .13;
    yfac = .13;
    set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+1.1 yfac*5]);
    set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+1.1 yfac*5]);
    
    for i = 1:2
        subplot(1,2,i)
        plot(ParamInVecMat(:,i)',ParamOutVecMat(:,i)','o'); hold on
        plot(0:.1:max(ParamInVecMat(:,i)),0:.1:max(ParamInVecMat(:,i)),':');
    end
      
end


