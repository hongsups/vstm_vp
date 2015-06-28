function create_modelfit_delay
%% run subject data
SbjNameCell = {'HJK','HS','JYP','RRS','YS'};
DelayIdxVec = [1 2 3 6];
Nvec = [1 2 4 6];
nParam = 3;
nSbj = length(SbjNameCell);
nDelay = length(DelayIdxVec);

if exist('de_delay/de_delay_vp_HJK_1000.mat','file')
    fprintf('Skipping model-fitting... \n');
    
else
    TestIdx = 0;
    for i = 1:nSbj
        for j = 1:nDelay
            SbjIdx = i;
            DelayIdx = DelayIdxVec(j);
            fit_de_delay_vp(SbjIdx,DelayIdx,TestIdx);
        end
    end
end

%% load all delay datas per subject and take average over subjects 
LLHMat = zeros(nDelay,nSbj);
ParamOutMat = zeros(nDelay,nSbj,nParam);
CSDMat = zeros(nSbj,length(Nvec),nDelay);
for i = 1:nDelay
    LLHVec = [];
    ParamOutVecMat = [];
    for j = 1:nSbj
        Delay = DelayIdxVec(i)*1000;
        SbjIdx = SbjNameCell{j};
        
        FileName = ['de_delay/de_delay_vp_' SbjIdx '_' num2str(Delay) '.mat'];
        load(FileName,'fval','ParamOutVec','SbjData');
        
        LLHVec = [LLHVec; fval];
        ParamOutVecMat = [ParamOutVecMat; ParamOutVec];
        for k = 1:length(Nvec)
            data_sbj = SbjData{Nvec(k),DelayIdxVec(i)};
            data_sbj = data_sbj(:,3)*pi/180;
            CSDMat(j,k,i) = compute_csd_de_ori(data_sbj);
        end
    end
    LLHMat(i,:) = LLHVec;
    ParamOutMat(i,:,:) = ParamOutVecMat;
end

LLHMat_mix = zeros(nDelay,nSbj);
ParamOutMat_mix = zeros(nDelay,nSbj,nParam-1);
for i = 1:nDelay
    LLHVec = [];
    ParamOutVecMat = [];
    for j = 1:nSbj
        Delay = DelayIdxVec(i)*1000;
        SbjIdx = SbjNameCell{j};
        
        FileName = ['de_delay/de_delay_mix_' SbjIdx '_' num2str(Delay) '.mat'];
        load(FileName,'fval','ParamOutVec');
        
        LLHVec = [LLHVec; fval];
        ParamOutVecMat = [ParamOutVecMat; ParamOutVec];
    end
    LLHMat_mix(i,:) = LLHVec;
    ParamOutMat_mix(i,:,:) = ParamOutVecMat;
end

figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = .13;
yfac = .13;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+1 yfac*5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+1 yfac*5]);

for i = 1:nParam-1
    ParamMat = squeeze(ParamOutMat_mix(:,:,i));
    
    subplot(1,nParam-1,i)
    plot(1:nDelay,mean(ParamMat,2),':o');
    hold all
    h = errorbar(1:nDelay,mean(ParamMat,2),std(ParamMat,0,2)/sqrt(nSbj));
    set(h,'linestyle','none');
    if i == 1
        title('J');
    else
        title('w');
    end
    set(gca,'Xtick',1:nDelay);
    set(gca,'XTickLabel',{'1','2','3','6'});
    xlim([0.5 4.5]);
    
    if i == 2
        ylim([0 1]);
    elseif i == 1
        ylim([0 50]);
    end
    
%     anova2(ParamMat)
    
end

%% plot CSD
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = .13;
yfac = .13;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+.5 yfac*5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+.5 yfac*5]);
meanCSD = squeeze(mean(CSDMat,1));
semCSD = squeeze(std(CSDMat,0,1))/sqrt(length(SbjNameCell));
plot(meanCSD',':o'); hold on
errorbar(meanCSD',semCSD');
title('CSD');
legend({'N=1','N=2','N=4','N=6'},'Location','Best');
xlabel('Delay (s)');
ylabel('CSD (rad)');
xlim([0.5 4.5]);
ylim([0 pi/4]);
set(gca,'XTick',1:4);set(gca,'XTickLabel',[1 2 3 6]);
set(gca,'YTick',0:pi/8:pi/4);set(gca,'YTickLabel');

%% CSD anova: 5 (sub) 4 (Nvec) 4 (delay)
CSDMat_for_anova = reshape(CSDMat,[],length(Nvec)); % 5 4 4 
reps = length(SbjNameCell);
p = anova2(CSDMat_for_anova,reps)

%% plot LLH, Jbar, power, and tau
% scatter plot with mean and sem among subjects
% plot MLE of each parameter against different delays

%% LLH
% figure;
% set(gca,'FontSize',11,'FontName','Arial');
% xfac = .13;
% yfac = .13;
% set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+.5 yfac*5]);
% set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+.5 yfac*5]);
% 
% bar(mean(LLHMat,2),'FaceColor','y');
% xlabel('Delays (s)');ylabel('Model likelihood');
% hold all
% h = errorbar(1:nDelay,mean(LLHMat,2),std(LLHMat')/sqrt(nSbj));
% set(h,'linestyle','none');
% set(gca,'Xtick',1:nDelay);
% set(gca,'XTickLabel',{'1','2','3','6'});

%% MLE
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = .13;
yfac = .13;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+2.1 yfac*5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+2.1 yfac*5]);

ParamNameCell = {'J_{bar}','{\alpha}','{\tau}','CV'};

for i = 1:length(ParamNameCell)
    
    if i < 4
        ParamMat = squeeze(ParamOutMat(:,:,i));
    else
        ParamMat = squeeze(sqrt(ParamOutMat(:,:,3)./ParamOutMat(:,:,1)));
    end
    ParamName = ParamNameCell{i};
    
    subplot(1,length(ParamNameCell),i)
    plot(1:nDelay,mean(ParamMat,2),':o');
    hold all
    h = errorbar(1:nDelay,mean(ParamMat,2),std(ParamMat,0,2)/sqrt(nSbj));
    set(h,'linestyle','none');
    title(ParamName);
    set(gca,'Xtick',1:nDelay);
    set(gca,'XTickLabel',{'1','2','3','6'});
    xlim([0.5 4.5]);
    
    if i == 2
        ylim([0 1.5]);
    elseif i == 3
        ylim([0 20]);
    elseif i == 1
        ylim([0 50]);
    elseif i == 4
        ylim([0.5 1.5]);
        text(.5,8,num2str(mean(ParamMat,2)',2));
    end
end

% anova
for i = 1:length(ParamNameCell)
    if i < 4
        ParamMat = squeeze(ParamOutMat(:,:,i));
    else
        ParamMat = squeeze(sqrt(ParamOutMat(:,:,3)./ParamOutMat(:,:,1)));
    end

    anova2(ParamMat)
end
% distribution
colorvec = get(gca, 'ColorOrder');

figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = .13;
yfac = .13;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+1 yfac*15]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+1 yfac*15]);

meanJbar = mean(squeeze(ParamOutMat(:,:,1)),2);
meanTau = mean(squeeze(ParamOutMat(:,:,3)),2);

nsteps = 1000;
bins = linspace(0,100,nsteps);

gampdfmat = zeros(length(SbjNameCell),length(DelayIdxVec),length(bins));
for i = 1:length(SbjNameCell)
    subplot(length(SbjNameCell),1,i)
    for j = 1:length(DelayIdxVec)
       temp = gampdf(bins,ParamOutMat(j,i,1)/ParamOutMat(j,i,3),ParamOutMat(j,i,3));
       gampdfmat(i,j,:) = temp;
       plot(bins,temp,'Color',colorvec(j,:)); hold on
    end
    if i == 1
        legend({'1s','2s','3s','6s'},'Location','Best');
        title('Distribution of precision, J');
    end
    if i == length(SbjNameCell)
        xlabel('J');
    end
    ylabel('Probability');
end

figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = .13;
yfac = .13;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+1 yfac*5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+1 yfac*5]);

gampdf_sum = squeeze(sum(gampdfmat,1));
plot(bins,gampdf_sum);
legend({'1s','2s','3s','6s'},'Location','Best');
title('Distribution of precision J (all subjects)');
ylabel('Probability');
xlabel('J');

function csd2 = compute_csd_de_ori(radvec)

% 0. scael to [-pi, pi]
radvec2 = radvec*2;

% 1. compute the length
x = mean(sin(radvec2));
y = mean(cos(radvec2));
[~,R] = cart2pol(x,y);
    
% 2. csd
csd = sqrt(-2*log(R));

% 3. scale back to [-pi/2, pi/2]
csd2 = csd/2;