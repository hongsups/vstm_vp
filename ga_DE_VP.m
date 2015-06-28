function ga_DE_VP(feature_idx,rel_idx,test_idx)








% this is for mturk








%%%%%%%%%%%%%%%%%%%%% 200??????????????????? WHY??????????



if ~exist('test_idx','var');
    test_idx = 0;
end

% feature_idx: 1-ori, 2-col
% rel_idx: 1-rel, 2-irr
%     1,1 -> 2
%     1,2 -> 1
%     2,1 -> 2
%     2,2 -> 1
features = {'ori', 'col'};
relevance = {'rel', 'irr'};

load_idx = mod(10*feature_idx + rel_idx,10)-1;
if load_idx == 0
    load_idx = 2;
end

% load subject data here and use it as a passing input for ga
if test_idx == 1 % fake data test
    fname = ['ga_DE_VP_test_' features{feature_idx} '_' relevance{rel_idx} '_' num2str(sbj_idx) '.mat'];
    load(fname);
elseif test_idx == 0
    % load subject data
    % load data
    [error_cell{1}, error_cell{2}] = loaddata_de(load_idx);
    error_vec_input = error_cell{feature_idx};
    data_sbj = error_vec_input*pi/180;
    if feature_idx == 2 % color
        data_sbj = data_sbj/2;
    end
    error_vec_input = abs(data_sbj);
    fname = ['ga_DE_VP' features{feature_idx} '_' relevance{rel_idx}];
end

% set population size (default: 20, but with small # of params, 10 is enough)
% set initial range (default is [0 1]);
% Handle to the function that produces mutation children
% Mutation function: default (Gaussian)
% StallGenLimit = ## (after ## gen of stall, it stops)

popSize = 64;
stallGenLim = 20;
genNum = 64;

LB = [0 0];
UB = [Inf Inf];
initRange = [1; 60];
nvar = 2;
    
% plot the process
opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping});
% set popSize, initial range, stallgenlim, gennum
opts = gaoptimset(opts,'PopulationSize',popSize,'PopInitRange',initRange,'StallGenLimit',stallGenLim,'Generations',genNum);

% run ga solver: ub, lb, passing subject index as an additional variable
FitnessFun = {@LLH,error_vec_input};
[x, fval, ~, Output] = ga(FitnessFun,nvar,[],[],[],[],LB,UB,[],[],opts);

fprintf('The number of generations was : %d\n', Output.generations);
fprintf('The number of function evaluations was : %d\n', Output.funccount);
fprintf('The best function value found was : %g\n', -fval);
fprintf('The MLE of Jbar is : %g\n', x(1));
fprintf('The MLE of tau is : %g\n', x(2));

% LLH input: error_vec_ori_irr,error_vec_col_rel,error_vec_ori_rel,error_vec_col_irr

% MLE (output)
Jbar_out = x(1);
tau_out = x(2);

param_out_vec = x;
if test_idx == 1
    param_in_vec = [Jbar_in tau_in];
elseif test_idx == 0
    param_in_vec = [];
end

% generate fake data
ntrials = 15000;
data_model = create_fakedata_DE(ntrials, Jbar_out, tau_out);

% model fitting
generate_modelfit_DE_VP(fname,data_sbj,data_model,param_in_vec,param_out_vec,fval);

function LLH_ga_VP = LLH(param_vec,error_vec_input)

Jbar_out = param_vec(1);
tau_out = param_vec(2);

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

% MC simulation
sampleSize = 1280;
error_vec = linspace(0, pi/2, 200);

J = gamrnd(Jbar_out/tau_out, tau_out, 1, sampleSize);
kappa = interp1(J_vec,k_vec,J);
kappa = min(kappa,700);
denom = 2*pi*besseli(0,kappa);

% generating the mixture of VM
p_error = zeros(1,length(error_vec));
for ie = 1:length(error_vec)
    err = error_vec(ie);

    p_error(ie) = mean(1./denom.*exp(kappa.*cos(err)));
end

error_idx = zeros(1,length(error_vec));
for i = 1:length(error_vec_input)
    error_idx(i) = interp1(error_vec,1:length(error_vec),error_vec_input(i),'nearest','extrap'); % error w.r.t. target
end

% compute likelihood: get p(response | data) for all trials
p_resp = p_error(error_idx); 

% avoid log(0)
p_resp = max(p_resp,eps); 

p_LLH = sum(log(p_resp));

% since GA FINDS MINIMUM, we need to flipt the sign.
LLH_ga_VP = -p_LLH;

function generate_modelfit_DE_VP(fname,data_sbj,data_model,param_in_vec,param_out_vec,fval)

% ntrials_rel = 600;
% ntrials_irr = 18000;
nbins = 1000;

mean_pdf_error_model = mean(data_model,2);

bins = linspace(-pi/2,pi/2,nbins);

nbins_hist = 50;
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = .13;
yfac = .13;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+.5 yfac*5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+.5 yfac*5]);

[f1,x1] = hist(data_sbj,nbins_hist);
bar(x1,f1/trapz(x1,f1));
title(fname); hold on
plot(bins, mean_pdf_error_model,'r');

text(.8,2.2,['J_{bar} = ' num2str(param_out_vec(1),4)]);
text(.8,2,['tau = ' num2str(param_out_vec(2),4)]);

set(gca,'XTick',-pi/2:pi/2:pi/2);
set(gca,'XTickLabel',[-90;0;90]);
xlabel('Estimation error');
xlim([-2/pi-1.5 2/pi+1.5]);
ylim([0 3]);
% print(gcf,'-dpng',[fname '_ga' '.png']);
saveas(gcf,[fname '_ga'],'fig');

fname2 = [fname '.mat'];
save(fname2);

% 4. draw figure: YOU NEED TO NORMALIZE HISTOGRAM BY SUM

% nbin = 50;
% [f,x] = hist(error_vec,nbin);
% n_cnt = histc(error_vec,x);
% plot(x,n_cnt);

% bar(x,f/trapz(x,f));hold on
% plot(x,g,'r');hold off




