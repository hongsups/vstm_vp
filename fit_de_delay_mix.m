function fit_de_delay_mix(SbjIdx,DelayIdx,TestIdx)

% providing summary statistics for DE: w and CSD
% assuming that error distribution is mixture of uniform (1-w) and VM (w)
% parameters: J, power, w

if ~exist('DelayIdx','var');
    DelayIdx = 0;
    TestIdx = 1;
end

% set size (N): 1,2,4,6
% delay (T): 1000,2000,3000,6000 ms
% First, we fit data nonparametrically in terms of T

% load data
PathName = 'de_delay/';
if TestIdx == 1 % fake data test
%     FileName = [PathName 'de_delay_mix_test_' num2str(SbjIdx) '.mat'];
    FileName = [PathName 'de_delay_mix_test_' num2str(SbjIdx)];
    SbjName = ['T' num2str(SbjIdx)];
    load([FileName '.mat']);
elseif TestIdx == 0
    SbjNameCell = {'HJK','HS','JYP','RRS','YS'};
    SbjName = SbjNameCell{SbjIdx};
%     FileName = [PathName 'de_delay_mix_' SbjName '_' num2str(DelayIdx*1000) '.mat'];
    FileName = [PathName 'de_delay_mix_' SbjName '_' num2str(DelayIdx*1000)];
    SbjData = load_de_delay(SbjName);
end

% GA settings
% popSize = 10;
% genNum = 10;

popSize = 50;
genNum = 30;
stallGenLim = 20;

LB = [0 0];
UB = [100 1];
initRange = [LB; UB];
nvar = 2;   % J, power, w
    
% plot the process
opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping});
% set popSize, initial range, stallgenlim, gennum
opts = gaoptimset(opts,'PopulationSize',popSize,'PopInitRange',initRange,'StallGenLimit',stallGenLim,'Generations',genNum);

% run ga solver: ub, lb, passing subject index as an additional variable
% FitnessFun = {@LLH,error_vec_input};
FitnessFun = {@LLH,TestIdx,SbjData,DelayIdx};
[x, fval, ~, Output] = ga(FitnessFun,nvar,[],[],[],[],LB,UB,[],[],opts);

fprintf('The number of generations was : %d\n', Output.generations);
fprintf('The number of function evaluations was : %d\n', Output.funccount);
fprintf('The best function value found was : %g\n', -fval);
fprintf('The MLE of J is : %g\n', x(1));
fprintf('The MLE of w is : %g\n', x(2));

% LLH input: error_vec_ori_irr,error_vec_col_rel,error_vec_ori_rel,error_vec_col_irr

% MLE (output)
JOut = x(1);
wOut = x(2);

ParamOutVec = x;
if TestIdx == 0
    ParamInVec = [];
end

% generate fake data
Nvec = [1 2 4 6];

DataModelCell = cell(1,4);
for i = 1:length(Nvec)
    DataModelCell{i} = create_fakedata_DE_mix(JOut, 1, wOut, Nvec(i));
end

% model fitting
generate_modelfit_DE_mix(TestIdx,SbjName,FileName,SbjData,DelayIdx,DataModelCell,ParamInVec,ParamOutVec,fval);

function LLH_ga_mix = LLH(ParamVec,TestIdx,SbjData,DelayIdx)

JOut = ParamVec(1);
PowerOut = 1;
wOut = ParamVec(2);

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

Nvec = [1 2 4 6];

p_LLH_temp = zeros(1,length(Nvec));
for i = 1:length(Nvec)

    N = Nvec(i);
    
    % load subject data
    if TestIdx == 0
        SbjDataMat = SbjData{N,DelayIdx};
%         ErrorVecInput = SbjDataMat(:,3)*pi/180*2;
%         ErrorVecInput = abs(SbjDataMat(:,3))*pi/180;
        ErrorVecInput = SbjDataMat(:,3)*pi/180;
    else
        ErrorVecInput = SbjData{i};
%         ErrorVecInput = abs(SbjData{i});
    end

    % sort data based on set sizes
    J = JOut*N^(-PowerOut);
    kappa = interp1(J_vec,k_vec,J);
    kappa = min(kappa,700);
    denom = 2*pi*besseli(0,kappa);

    ErrorVec = linspace(-pi/2, pi/2, length(ErrorVecInput));
%     ErrorVec = linspace(0, pi, length(ErrorVecInput));
    p_error = 1./denom.*exp(kappa.*cos(ErrorVec))*wOut + (1-wOut)/pi;

    error_idx = zeros(1,length(ErrorVec));
    for ik = 1:length(ErrorVecInput)
        error_idx(ik) = interp1(ErrorVec,1:length(ErrorVec),ErrorVecInput(ik),'nearest','extrap'); % error w.r.t. target
    end

    % compute likelihood: get p(response | data) for all trials
    p_resp = p_error(error_idx); 

    % avoid log(0)
    p_resp = max(p_resp,eps); 

    p_LLH_temp(i) = sum(log(p_resp));
end
    
% p_LLH = sum(log(p_resp));
p_LLH = sum(p_LLH_temp);

% since GA FINDS MINIMUM, we need to flipt the sign.
LLH_ga_mix = -p_LLH;

function generate_modelfit_DE_mix(TestIdx,SbjName,FileName,SbjData,DelayIdx,DataModelCell,ParamInVec,ParamOutVec,fval);

nbins = 1000;
Nvec = [1 2 4 6];

figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = .13;
yfac = .13;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+2 yfac*5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+2 yfac*5]);

% separate model figure
for i = 1:length(Nvec)
    N = Nvec(i);
    
    if TestIdx == 0
        data_sbj = SbjData{N,DelayIdx};
        data_sbj = data_sbj(:,3)*pi/180*2;
        data_model = DataModelCell{i};
    else
        FileName = ['de_delay/' SbjName '-delay-mix-errorfit'];
        data_sbj = SbjData{i};
        data_model = DataModelCell{i};
    end        
    
    mean_pdf_error_model = mean(data_model,2);

    bins = linspace(-pi,pi,nbins);
    nbins_hist = 12;
    subplot(1,length(Nvec),i)
    
    [f1,x1] = hist(data_sbj,nbins_hist);
    bar(x1,f1/trapz(x1,f1));
    title([SbjName ' ' num2str(DelayIdx*1000) ' ' 'N=' num2str(N)]); hold on
    
    plot(bins, mean_pdf_error_model,'r');

    CSD = compute_csd_de_ori(data_sbj); % unit: radian
    
    [f1,x1] = hist(data_sbj,nbins_hist);
    bar(x1,f1/trapz(x1,f1));
    title([SbjName ' ' num2str(DelayIdx*1000) ' ' 'N=' num2str(N)]); hold on
    
    plot(bins, mean_pdf_error_model,'r');

    xcord = -pi*4/3;
    ycord = 2.0;
    if i == length(Nvec)
        text(xcord,ycord+0.2,['LLH:' num2str(-fval,4)]);
        text(xcord,ycord,['Out:' num2str(ParamOutVec,4)]);
        if TestIdx == 1
            text(xcord,ycord-0.2,['In:' num2str(ParamInVec,4)]);
        end
    end
    
    set(gca,'XTick',-pi:pi:pi);
    set(gca,'XTickLabel',[-90;0;90]);
    xlabel('Estimation error');
    xlim([-pi-1.5 pi+1.5]);
    ylim([0 3]);
    text(xcord,ycord-0.4,num2str(CSD,4));
    
    % print(gcf,'-dpng',[FileName '_ga' '.png']);
end

saveas(gcf,[FileName '.fig'],'fig');
fname2 = [FileName '.mat'];
save(fname2);

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
