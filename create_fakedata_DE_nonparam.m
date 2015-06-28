function [error_pdf_mat] = create_fakedata_DE_nonparam(ntrials, J_bar, tau)

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

%% method 1: pdf based on error histogram (large sample size)
% error_vec = zeros(ntrials,1);
% trialcnt = 0;
% while (trialcnt<ntrials)
%     trialcnt = trialcnt + 1;
%     
%     J = gamrnd(J_bar/tau, tau);
%     kappa = interp1(J_vec,k_vec,J);
%     kappa = min(kappa,700);
% 
%     error_vec(trialcnt) = circ_vmrnd(0,kappa);
% end
%
% histogram
% nbin = 50;
% [f,x] = hist(error_vec,nbin);
% n_cnt = histc(error_vec,x);
% plot(x,n_cnt);

%% method 2: averaging distributions from Js
Js = gamrnd(J_bar/tau, tau, 1, ntrials);
kappas = interp1(J_vec,k_vec,Js);
kappas = min(kappas,700);
nbin = 1000;
bins = linspace(-pi,pi,nbin);
error_pdf_mat = zeros(nbin, ntrials);
for i = 1:ntrials
    k = kappas(i);
    [error_pdf_mat(:,i),~] = circ_vmpdf(bins,0,k);
end
% plot pdf of the mixture
% figure;
% set(gca,'FontSize',11,'FontName','Arial');
% xfac = .70;
% yfac = .13;
% set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac+.5 yfac*5]);
% set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac+.5 yfac*5]);
% 
% plot(mean(p,2));
