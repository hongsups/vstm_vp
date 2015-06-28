function csd_de

WRONG CSD

% subject names (N = 5)
SbjNameCell = {'HJK','HS','JYP','RRS','YS'};
Tvec = [1 2 3 6];
Nvec = [1 2 4 6];

csdMat = zeros(length(SbjNameCell),length(Tvec),length(Nvec));
for i = 1:length(SbjNameCell)
    
    % load subject data: /de_delay/de_delay_vp_SbjName_1000.mat
    % SbjData: 6x6 cell (rows: set sizes (1,2,4,6), columns: delays (1,2,3,6))
    % SbjData(i,j) will load 120x3 matrix
    % Here, the 3rd column has the error (=response - stimulus)
    SbjName = SbjNameCell{i};
    load(['de_delay/de_delay_vp_' SbjName '_1000.mat'],'SbjData');
    
    for j = 1:length(Nvec)
        N = Nvec(j);
        for k = 1:length(Tvec)
            T = Tvec(k);
            degvec = SbjData{N,T}(:,3);
            
            csdMat(i,j,k) = compute_csd_de_ori(degvec);
        end
    end    
end

% plot: average
% x: T, y: csd, w/ different Ns
mean_csd = squeeze(mean(csdMat,1)); % Dim1: N, Dim2: T
sem_csd = squeeze(std(csdMat))/sqrt(length(SbjNameCell));

save all
qdas

figure;
set(gca,'FontSize',11,'FontName','Arial'); xfac = .05; yfac = .23;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+.5 yfac*2.5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+.5 yfac*2.5]);

for i = 1:length(Nvec)
    plot(Tvec,mean_csd(i,:),'o'); hold on
    errorbar(Tvec, mean_csd(i,:), sem_csd(i,:));
    ylim([0 pi/4]);%set(gca,'YTick',20:20:60);
    xlim([.5 4.5]);set(gca,'XTick',1:4);
    xlabel('T (sec)'); 
    if i == 1
        ylabel('CSD (deg)');
    end
    hold on
end
% savefig('csd_de_all.fig');
print(gcf,'-dpng','csd_de_all.png');

% 2-way anova: row factor / column factor / and replicate (dimensions)
% change the dimension order
csdMat_for_anova = reshape(csdMat,[],length(Nvec));
reps = length(SbjNameCell);
p = anova2(csdMat_for_anova,reps)

% % plot: individual
% for i = 1:length(SbjNameCell)
%     figure;
%     set(gca,'FontSize',11,'FontName','Arial'); xfac = .05; yfac = .23;
%     set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+.5 yfac*2.5]);
%     set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+.5 yfac*2.5]);
%     SbjName = SbjNameCell{i};
%     for j = 1:length(Nvec)
%         %subplot(1,length(Nvec),j);
%         plot(Tvec,squeeze(csdMat(i,j,:)),':o');
%         ylim([0 pi/2]);%set(gca,'YTick',20:20:60);
%         xlim([.5 4.5]);set(gca,'XTick',1:4);
%         title([SbjName ' N=' num2str(Nvec(j))]);
%         xlabel('T (sec)'); 
%         if j == 1
%             ylabel('CSD (deg)');
%         end
%         hold on
%     end
% %     savefig(['csd_de_' SbjName '.fig']);
%     print(gcf,'-dpng',['csd_de_' SbjName '.png']);
% end

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
