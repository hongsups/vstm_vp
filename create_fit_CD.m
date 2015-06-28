function [hit_rate false_alarm_rate prop_magN m_hit_rate m_false_alarm_rate m_prop_magN] = create_fit_CD(data_sbj,data_model,expid)

nBins = 9;
bin = linspace(0,90,nBins+1);
bin = bin(1:end-1)+diff(bin(1:2))/2;
binsize = diff(bin(1:2));

if expid == 1
    Nvec_exp = 2:2:8;                       % set sizes used in the experiments
else
    Nvec_exp = [1 2 4 8];
end

% subject data
response_idx_vec = data_sbj(:,2);
N_idx_vec = data_sbj(:,5);
if max(data_sbj(:,1))==8        % real subject data
    change_idx_vec = ceil((data_sbj(:,1))/max(data_sbj(:,1)));
    mag_vec = data_sbj(:,6)/pi*180.*change_idx_vec;     % [0,180] (You need to multiply change_idx)
else
    change_idx_vec = data_sbj(:,1);
    mag_vec = abs(data_sbj(:,6));
end
mag_vec(mag_vec>90)=180-mag_vec(mag_vec>90);        % remapping the delta vector to [0,90]

% model data
m_response_idx_vec = data_model(:,2);
m_N_idx_vec = data_model(:,5);
if max(data_model(:,1))==8        % real subject data
    m_change_idx_vec = ceil((data_model(:,1))/max(data_model(:,1)));
    m_mag_vec = data_model(:,6)/pi*180.*m_change_idx_vec;     % [0,180] (You need to multiply change_idx)
else
    m_change_idx_vec = data_model(:,1);
    m_mag_vec = abs(data_model(:,6));
end
m_mag_vec(m_mag_vec>90)=180-m_mag_vec(m_mag_vec>90);    % remapping the delta vector to [0,90]

% hit, false alarm
hit_rate = zeros(1,length(Nvec_exp));
false_alarm_rate = zeros(1,length(Nvec_exp));
prop_magN = zeros(length(Nvec_exp),length(bin));
m_hit_rate = zeros(1,length(Nvec_exp));
m_false_alarm_rate = zeros(1,length(Nvec_exp));
m_prop_magN = zeros(length(Nvec_exp),length(bin));
for i = 1:length(Nvec_exp)
    N = Nvec_exp(i);
    hit_rate(i) = sum((response_idx_vec == 1) & (change_idx_vec == 1) & (N_idx_vec == N))/sum(change_idx_vec == 1 & N_idx_vec == N);
    false_alarm_rate(i) = sum((response_idx_vec == 1) & (change_idx_vec == 0) & (N_idx_vec == N))/sum(change_idx_vec == 0 & N_idx_vec == N);

    m_hit_rate(i) = sum((m_response_idx_vec == 1) & (m_change_idx_vec == 1) & (m_N_idx_vec == N))/sum(m_change_idx_vec == 1 & m_N_idx_vec == N);
    m_false_alarm_rate(i) = sum((m_response_idx_vec == 1) & (m_change_idx_vec == 0) & (m_N_idx_vec == N))/sum(m_change_idx_vec == 0 & m_N_idx_vec == N);

    for jj = 1:length(bin)
        ind_sbj = (abs(mag_vec) >= bin(jj)-binsize/2) & (abs(mag_vec) < bin(jj)+binsize/2) & N_idx_vec == N & response_idx_vec==1 & change_idx_vec == 1;
        prop_magN(i,jj) = sum(ind_sbj)/sum(N_idx_vec == N & abs(mag_vec) >= bin(jj)-binsize/2 & abs(mag_vec) < bin(jj)+binsize/2 & change_idx_vec == 1);
        
        m_ind_sbj = (abs(m_mag_vec) >= bin(jj)-binsize/2) & (abs(m_mag_vec) < bin(jj)+binsize/2) & m_N_idx_vec == N & m_response_idx_vec==1 & m_change_idx_vec == 1;
        m_prop_magN(i,jj) = sum(m_ind_sbj)/sum(m_N_idx_vec == N & abs(m_mag_vec) >= bin(jj)-binsize/2 & abs(m_mag_vec) < bin(jj)+binsize/2 & m_change_idx_vec == 1);
    end    
end

figure;
set(gca,'FontSize',11,'FontName','Arial');
fac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 2 2]*fac);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 2 2]*fac);
% performance by set size
subplot(1,2,1); 
plot(Nvec_exp, hit_rate,'ob');
hold on;
plot(Nvec_exp, false_alarm_rate,'or');
hold on;
plot(Nvec_exp, m_hit_rate,':b');
hold on;
plot(Nvec_exp, m_false_alarm_rate,':r');
% plot(Nvec_exp, m_prop_response_N,':'); 
ylabel('Proportion reports "change"');axis([2 8 0 1]);
legend('hit rate','false alarm rate',1);

% performance by change magnitude: add 0 on x-axis for FA rate
bin_new = [0 bin];
prop_magN_new = [false_alarm_rate' prop_magN];
m_prop_magN_new = [m_false_alarm_rate' m_prop_magN];

subplot(1,2,2); 
plot(bin_new, prop_magN_new, 'o');
hold on
plot(bin_new, m_prop_magN_new, ':'); ylabel('Proportion correct');
axis([0 90 0 1]);legend(strcat('N= ',int2str(Nvec_exp')), 'Location','Best');
