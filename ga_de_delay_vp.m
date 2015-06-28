function ga_de_delay_vp(SbjIdx,DelayIdx,NIdx,TestIdx)

% for now just go to fit_de_delay_vp.m: I think we fixed the bug there!
% Whee!

if ~exist('DelayIdx','var');
    DelayIdx = 0;
    TestIdx = 1;
end

% here we fit de-delay data to VP model nonparametrically for set size and
% delay

% set size (N): 1,2,4,6
% delay (T): 1000,2000,3000,6000 ms
% First, we fit data nonparametrically in terms of T

% load data
PathName = 'de_delay/';
if TestIdx == 1 % fake data test
%     FileName = [PathName 'de_delay_vp_test_' num2str(SbjIdx) '.mat'];
    FileName = [PathName 'de_delay_vp_test_' num2str(SbjIdx)];
    SbjName = ['T' num2str(SbjIdx)];
    load([FileName '.mat']);
elseif TestIdx == 0
    SbjNameCell = {'HJK','HS','JYP','RRS','YS'};
    SbjName = SbjNameCell{SbjIdx};
%     FileName = [PathName 'de_delay_vp_' SbjName '_' num2str(DelayIdx*1000) '.mat'];
    FileName = [PathName 'de_delay_vp_' SbjName '_' num2str(DelayIdx*1000)];
    SbjDataAll = load_de_delay(SbjName);
    SbjData = SbjDataAll{NIdx,DelayIdx};
end

% GA settings
popSize = 10;
genNum = 10;
SampleSize = 50;

% popSize = 50;
% genNum = 30;
% SampleSize = 500;
stallGenLim = 20;

LB = [0 0];
UB = [100 100];
initRange = [LB; UB];
nvar = 2;
    
% plot the process
opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping});
% set popSize, initial range, stallgenlim, gennum
opts = gaoptimset(opts,'PopulationSize',popSize,'PopInitRange',initRange,'StallGenLimit',stallGenLim,'Generations',genNum);

% run ga solver: ub, lb, passing subject index as an additional variable
% FitnessFun = {@LLH,error_vec_input};
FitnessFun = {@LLH,TestIdx,SbjData,DelayIdx,NIdx,SampleSize};
[x, fval, ~, Output] = ga(FitnessFun,nvar,[],[],[],[],LB,UB,[],[],opts);

fprintf('The number of generations was : %d\n', Output.generations);
fprintf('The number of function evaluations was : %d\n', Output.funccount);
fprintf('The best function value found was : %g\n', -fval);
fprintf('The MLE of Jbar is : %g\n', x(1));
fprintf('The MLE of tau is : %g\n', x(2));
% LLH input: error_vec_ori_irr,error_vec_col_rel,error_vec_ori_rel,error_vec_col_irr

% MLE (output)
JbarOut = x(1);
tauOut = x(2);

ParamOutVec = x;
if TestIdx == 0
    ParamInVec = [];
end

% generate fake data
ntrials = 120; % 1920 trials for all conditions (480 per condition)
DataModel = create_fakedata_DE_nonparam(ntrials, JbarOut, tauOut);

% model fitting
generate_modelfit_DE_VP(TestIdx,SbjName,FileName,SbjData,DelayIdx,NIdx,DataModel,ParamInVec,ParamOutVec,fval);

function LLH_ga_VP = LLH(ParamVec,TestIdx,SbjData,DelayIdx,NIdx,SampleSize)

JbarOut = ParamVec(1);
tauOut = ParamVec(3);
N = NIdx;

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);
    
% load subject data
if TestIdx == 0
    SbjDataMat = SbjData{N,DelayIdx};
    ErrorVecInput = abs(SbjDataMat(:,3))*pi/180;
else
    ErrorVecInput = abs(SbjData(:,3));
end

J = gamrnd(JbarOut/tauOut, tauOut, 1, SampleSize);
kappa = interp1(J_vec,k_vec,J);
kappa = min(kappa,700);
denom = 2*pi*besseli(0,kappa);

%     ErrorVec = linspace(0, pi, length(SbjDataMat));
ErrorVec = linspace(0, pi, length(ErrorVecInput));
p_error = zeros(1,length(ErrorVec));
for ie = 1:length(ErrorVec)
    err = ErrorVec(ie);
    p_error(ie) = mean(1./denom.*exp(kappa.*cos(err)));
end

error_idx = zeros(1,length(ErrorVec));
for ik = 1:length(ErrorVecInput)
    error_idx(ik) = interp1(ErrorVec,1:length(ErrorVec),ErrorVecInput(ik),'nearest','extrap'); % error w.r.t. target
end

% compute likelihood: get p(response | data) for all trials
p_resp = p_error(error_idx); 

% avoid log(0)
p_resp = max(p_resp,eps); 

p_LLH_temp(i) = sum(log(p_resp));
    
p_LLH = sum(log(p_resp));

% since GA FINDS MINIMUM, we need to flipt the sign.
LLH_ga_VP = -p_LLH;

function generate_modelfit_DE_VP(TestIdx,SbjName,FileName,SbjData,DelayIdx,NIdx,DataModelCell,ParamInVec,ParamOutVec,fval);

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
        data_sbj = data_sbj(:,3)*pi/180;
        data_model = DataModelCell{i};
    else
        FileName = ['de_delay/' SbjName '-delay-errorfit'];
        data_sbj = SbjData{i};
        data_model = DataModelCell{i};
    end        
    
    mean_pdf_error_model = mean(data_model,2);

    bins = linspace(-pi,pi,nbins);
    nbins_hist = 25;
    subplot(1,length(Nvec),i)
    
    [f1,x1] = hist(data_sbj,nbins_hist);
    bar(x1,f1/trapz(x1,f1));
    title([SbjName ' ' num2str(DelayIdx*1000) ' ' 'N=' num2str(N)]); hold on
    
    plot(bins, mean_pdf_error_model,'r');

    if i == length(Nvec)
        xcord = -pi*4/3;
        ycord = 2.0;
        text(xcord,ycord+0.2,['LLH:' num2str(-fval,4)]);
        text(xcord,ycord,['Out:' num2str(ParamOutVec,4)]);

        if TestIdx == 1
            text(xcord,ycord-0.2,['In:' num2str(ParamInVec,4)]);
        end
    end
    
    set(gca,'XTick',-pi:pi:pi);
    set(gca,'XTickLabel',[-180;0;180]);
    xlabel('Estimation error');
    xlim([-pi-1.5 pi+1.5]);
    ylim([0 3]);
    % print(gcf,'-dpng',[FileName '_ga' '.png']);
end

saveas(gcf,[FileName '.fig'],'fig');
fname2 = [FileName '.mat'];
save(fname2);
