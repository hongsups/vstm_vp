function compareLLH_de_delay_np

% de delay data
% compare model VP and mixture model

%% load files
Delayvec = [1 2 3 6];
Nvec = [1 2 4 6];
PathName = 'de_delay/';
SbjNameCell = {'HJK','HS','JYP','RRS','YS'};
nvar = 2;

ParamOutVecMat_VP = zeros(length(Delayvec),length(Nvec),length(SbjNameCell),nvar);
LLHMat_VP = zeros(length(Delayvec),length(Nvec),length(SbjNameCell));
ParamOutVecMat_mix = zeros(length(Delayvec),length(Nvec),length(SbjNameCell),nvar);
LLHMat_mix = zeros(length(Delayvec),length(Nvec),length(SbjNameCell));

for i = 1:length(Delayvec)
    DelayIdx = Delayvec(i);
    for j = 1:length(Nvec)
        NIdx = Nvec(j);
        for k = 1:length(SbjNameCell)
            SbjName = SbjNameCell{k};
            FileName_VP = [PathName 'de_delay_vp_np_' SbjName '_' num2str(DelayIdx*1000) '_' num2str(NIdx) '.mat'];

            % load VP
            load(FileName_VP,'ParamOutVec','fval');
            ParamOutVecMat_VP(i,j,k,:) = ParamOutVec;
            LLHMat_VP(i,j,k) = fval;
            clear ParamOutVec; clear fval;

            % load mix
            FileName_mix = [PathName 'de_delay_mix_np_' SbjName '_' num2str(DelayIdx*1000) '_' num2str(NIdx) '.mat'];
            load(FileName_mix,'ParamOutVec','fval');
            ParamOutVecMat_mix(i,j,k,:) = ParamOutVec;
            LLHMat_mix(i,j,k) = fval;
        end
    end
end

AICMat_VP = 2*nvar+2*LLHMat_VP;
AICMat_mix = 2*nvar+2*LLHMat_mix;
% meanAIC_VP = mean(AICMat_VP,3);
% meanAIC_mix = mean(AICMat_mix,3);
% semAIC_VP = std(AICMat_VP,0,3)/sqrt(length(SbjNameCell));
% semAIC_mix = std(AICMat_mix,0,3)/sqrt(length(SbjNameCell));

diffAICMat = AICMat_mix - AICMat_VP;
mean_diffAIC = mean(diffAICMat,3);
sem_diffAIC = std(diffAICMat,0,3)/sqrt(length(SbjNameCell));

figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = .13;
yfac = .13;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+2 yfac*5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+2 yfac*5]);

for i = 1:length(mean_diffAIC)
    subplot(1,4,i);
    Nvec = [1 2 4 6];
    N = Nvec(i);
    
    bar(mean_diffAIC(:,i),'FaceColor','y');
    ylim([-3 3]);xlim([0 5]);
    ylabel('Relative AIC to VP');
    title(['N=' num2str(N) ' AIC (=mixture-VP)']);

    hold all
    h = errorbar(1:length(mean_diffAIC),mean_diffAIC(:,i),sem_diffAIC(:,i));
    set(h,'linestyle','none');
    set(gca,'Xtick',1:length(mean_diffAIC));
    set(gca,'XTickLabel',{'1','2','3','6'});
    xlabel('Delay (sec)')
end
print(gcf,'-dpng','de_delay_np_relative_AIC.png');
% saveas(gcf,['BMC_CD_group_K_' exp_name '.emf'],'emf');

