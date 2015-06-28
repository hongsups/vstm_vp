function compareLLH_de_delay

% de delay data
% compare model VP and mixture model

%% load files
Delayvec = [1 2 3 6];
PathName = 'de_delay/';
SbjNameCell = {'HJK','HS','JYP','RRS','YS'};
nvarVP = 3;
nvarmix = 2;
nSbj = length(SbjNameCell);

ParamOutVecMat_VP = zeros(length(Delayvec),length(SbjNameCell),nvarVP);
LLHMat_VP = zeros(length(Delayvec),length(SbjNameCell));
ParamOutVecMat_mix = zeros(length(Delayvec),length(SbjNameCell),nvarmix);
LLHMat_mix = zeros(length(Delayvec),length(SbjNameCell));

for i = 1:length(Delayvec)
    DelayIdx = Delayvec(i);
    for j = 1:length(SbjNameCell)
        SbjName = SbjNameCell{j};
        FileName_VP = [PathName 'de_delay_vp_' SbjName '_' num2str(DelayIdx*1000)];
        
        % load VP
        load(FileName_VP,'ParamOutVec','fval');
        ParamOutVecMat_VP(i,j,:) = ParamOutVec;
        LLHMat_VP(i,j) = fval;
        clear ParamOutVec; clear fval;
        
        % load mix
        FileName_mix = [PathName 'de_delay_mix_' SbjName '_' num2str(DelayIdx*1000)];
        load(FileName_mix,'ParamOutVec','fval');
        ParamOutVecMat_mix(i,j,:) = ParamOutVec;
        LLHMat_mix(i,j) = fval;
        
    end
end

AICMat_VP = 2*nvarVP+2*LLHMat_VP;
AICMat_mix = 2*nvarmix+2*LLHMat_mix;
meanAIC_VP = mean(AICMat_VP,2);
meanAIC_mix = mean(AICMat_mix,2);
semAIC_VP = std(AICMat_VP,0,2)/sqrt(length(SbjNameCell));
semAIC_mix = std(AICMat_mix,0,2)/sqrt(length(SbjNameCell));

diffAICMat = AICMat_mix - AICMat_VP;
mean_diffAIC = mean(diffAICMat,2);
sem_diffAIC = std(diffAICMat,0,2)/sqrt(length(SbjNameCell));

% parameter estimates:Jbar,power,tau // J,power,w
meanParamVP = zeros(length(Delayvec),nvarVP);
meanParammix = zeros(length(Delayvec),nvarmix);
semParamVP = zeros(length(Delayvec),nvarVP);
semParammix = zeros(length(Delayvec),nvarmix);
for i = 1:nvarVP
    meanParamVP(:,i) = mean(squeeze(ParamOutVecMat_VP(:,:,i)),2);
    semParamVP(:,i) = std(squeeze(ParamOutVecMat_VP(:,:,i)),0,2)/sqrt(nSbj);
end

for i = 1:nvarmix
    meanParammix(:,i) = mean(squeeze(ParamOutVecMat_mix(:,:,i)),2);
    semParammix(:,i) = std(squeeze(ParamOutVecMat_mix(:,:,i)),0,2)/sqrt(nSbj);
end

% figure;
% set(gca,'FontSize',11,'FontName','Arial');
% xfac = .13;
% yfac = .13;
% set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+1.6 yfac*10]);
% set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+1.6 yfac*10]);
% titlecellVP = {'Jbar','power','tau'};
% titlecellmix = {'J','power','w'};
% for i = 1:nvar
%     subplot(2,nvar,i);
%     plot(1:length(Delayvec),meanParamVP(:,i),':'); hold on
%     errorbar(1:length(Delayvec),meanParamVP(:,i),semParamVP(:,i),'o');
%     set(gca,'Xtick',1:length(Delayvec));
%     set(gca,'XTickLabel',{'1','2','3','6'});
%     xlim([.5 4.5]);
%     if i ~= 2
%         ylim([0 60]);
%     else
%         ylim([0 1]);
%     end
%     title(titlecellVP{i});
%     xlabel('Delay (sec)')
% 
%     subplot(2,nvar,i+nvar);
%     plot(1:length(Delayvec),meanParammix(:,i),':'); hold on
%     errorbar(1:length(Delayvec),meanParammix(:,i),semParammix(:,i),'o');
%     set(gca,'Xtick',1:length(Delayvec));
%     set(gca,'XTickLabel',{'1','2','3','6'});
%     xlim([.5 4.5]);
%     if i+nvar == 4
%         ylim([0 60]);
%     elseif i+nvar == 5
%         ylim([0 1.5]);
%     else
%         ylim([.5 1]);
%     end
%     title(titlecellmix{i});
%     xlabel('Delay (sec)')
% end
% 
% print(gcf,'-dpng','de_delay_power_param.png');

figure;
set(gca,'FontSize',11,'FontName','Arial');
fac = .5;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 2 2]*fac);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);
bar(mean_diffAIC,'FaceColor','y');
ylabel('Relative AIC to VP');
title('AIC difference (=mixture-VP) for each delay condition');

hold all
h = errorbar(1:length(mean_diffAIC),mean_diffAIC,sem_diffAIC);
set(h,'linestyle','none');
set(gca,'Xtick',1:length(mean_diffAIC));
set(gca,'XTickLabel',{'1','2','3','6'});
xlabel('Delay (sec)')

print(gcf,'-dpng','de_delay_power_relative_AIC.png');
saveas(gcf,'de_delay_power_relative_AIC.emf','emf');

