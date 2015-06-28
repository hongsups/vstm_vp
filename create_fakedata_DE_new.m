function [error_pdf_mat] = create_fakedata_DE_new(ntrials, J_bar, tau, power, N)

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

% averaging distributions from Js
Js = gamrnd(J_bar*N^(-power)/tau, tau, 1, ntrials);
kappas = interp1(J_vec,k_vec,Js);
kappas = min(kappas,700);
nbin = 1000;
bins = linspace(-pi,pi,nbin);
error_pdf_mat = zeros(nbin, ntrials);
for i = 1:ntrials
    k = kappas(i);
    [error_pdf_mat(:,i),~] = circ_vmpdf(bins,0,k);
end
