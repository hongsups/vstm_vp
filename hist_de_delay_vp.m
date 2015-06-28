function hist_de_delay_vp(plot_idx,nbins_hist)

if ~exist('plot_idx','var')
    plot_idx = 1;
end

if ~exist('nbins_hist','var')
    nbins_hist = 25;
end

% project: de-delay (various set sizes and delays)
% reference code: fit_de_delay_vp (generate individual fit from ga)
% goal: generate average histogram among set sizes, delays, or each

% nbins_hist = 25;
PathName = 'de_delay/';
SbjNameCell = {'HJK','HS','JYP','RRS','YS'};

%% Load data
Nvec = [1 2 4 6];
Delayvec = [1 2 3 6];

% subject
nSbj = length(SbjNameCell);
SbjDataCell = cell(1,nSbj);
for i = 1:nSbj
    SbjName = SbjNameCell{i};
    SbjData = load_de_delay(SbjName); % 6 x 6 cell (row: set size, column: delay)
    SbjDataCell{i} = SbjData;
end

% model: model data is saved as histogram data (nbins x ntrials)
% DataModelCell: 1xNvec cell (each cell entry is nbins x ntrails matrix)
% DataModelCell is from each delay
ModelDataCell = cell(nSbj,length(Delayvec));
for j = 1:nSbj
    for i = 1:length(Delayvec)
        DelayIdx = Delayvec(i);
        ModelFileName = [PathName 'de_delay_vp_' SbjName '_' num2str(DelayIdx*1000)];
        load(ModelFileName,'DataModelCell');
        ModelDataCell{j,i} = DataModelCell;
    end
end

%% Calculate hist-counts using histcounts(X,edges)
% Subject and model data have different form.
% Subject data: it needs to be binned for histogram. nbins_hist = 25;
% Model data: finely binned and computed already: plot(bins,data_model)
% For subject data, we need to interpolate...
sbj_data_binned_mat_cell = cell(length(Delayvec),length(Nvec));
model_data_binned_mat_cell = cell(length(Delayvec),length(Nvec));
nbins_model = length(DataModelCell{1});
bin_model = linspace(-pi,pi,nbins_model);
edges = linspace(-pi,pi,nbins_hist);
for i = 1:length(Delayvec)
    DelayIdx = Delayvec(i);
    
    figure;
    set(gca,'FontSize',11,'FontName','Arial');
    xfac = .13;
    yfac = .13;
    set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+2 yfac*5]);
    set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+2 yfac*5]);
    
    for j = 1:length(Nvec)
        NIdx = Nvec(j);
        sbj_data_binned_mat = zeros(nSbj,nbins_hist);
        model_data_binned_mat = zeros(nSbj,nbins_hist);
        for k = 1:nSbj
            SbjIdx = k;
            
            % load subject data
            data_sbj = SbjDataCell{SbjIdx}{NIdx,DelayIdx}(:,3)*pi/180;
            
            % draw histogram and find the "edge" (x axis values)
            % [f1,edges] = hist(data_sbj,nbins_hist);
            f1 = hist(data_sbj,edges);
            % bar(x1,f1/trapz(x1,f1));
            sbj_data_binned = f1/trapz(edges,f1);  % normalizing
            
            % interpolate model value according to the edge
            mean_pdf_error_model = mean(ModelDataCell{SbjIdx,i}{j},2);
            model_data_binned = interp1(bin_model,mean_pdf_error_model,edges,'neareast','extrap');
            % hold on;plot(edge,C,'ro');
            
            % save
            sbj_data_binned_mat(k,:) = sbj_data_binned;
            model_data_binned_mat(k,:) = model_data_binned;
        end
        
        % mean and sem
        mean_sbj = mean(sbj_data_binned_mat);
        mean_model = mean(model_data_binned_mat);
        sem_sbj = std(sbj_data_binned_mat)/sqrt(nSbj);
        sem_model = std(model_data_binned_mat)/sqrt(nSbj);
        
        residual = model_data_binned_mat - sbj_data_binned_mat;
        mean_residual = mean(residual);
        sem_residual = std(residual)/sqrt(nSbj);
        
        % plot 
        subplot(1,length(Nvec),j);
        if plot_idx == 1    % raw data
            errorbar(edges,mean_sbj,sem_sbj,':o'); hold on;
            errorbar(edges,mean_model,sem_model,':o'); 
            xlim([-pi-1.5 pi+1.5]);
            ylim([0 3]);
        elseif plot_idx == 2    % residual
            %plot(edges,mean_residual);
            errorbar(edges,mean_residual,sem_residual);
            xlim([-pi-1.5 pi+1.5]);
            ylim([-.8 .8]);
            endx`
        set(gca,'XTick',-pi:pi:pi);
        set(gca,'XTickLabel',[-180;0;180]);
        xlabel('Estimation error');
        title(['T=' num2str(DelayIdx*1000) ' ' 'N=' num2str(NIdx)]);
        
        sbj_data_binned_mat_cell{i,j} = sbj_data_binned_mat;
        model_data_binned_mat_cell{i,j} = model_data_binned_mat;
    end
end
% X = randn(1000,1);
% edges = [-5 -4 -2 -1 -0.5 0 0.5 1 2 4 5];
% N = histcounts(X,edges)
% [N,edges] = histcounts(X, 'Normalization', 'probability')