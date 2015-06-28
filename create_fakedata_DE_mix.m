function error_pdf = create_fakedata_DE_mix(J, power, w, N)

% if ~exist('J','var')
%     J = 10;
%     power = 1;
%     w = 0.2;
%     N = 1;
% end

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

kappa = interp1(J_vec,k_vec,J*N^-(power));
kappa = min(kappa,700);
nbin = 1000;
bins = linspace(-pi,pi,nbin);

% error_pdf = circ_vmpdf(bins,0,kappa)*w + (1-w)/(2*pi);
error_pdf = circ_vmpdf(bins,0,kappa)*w + (1-w)/2/pi;