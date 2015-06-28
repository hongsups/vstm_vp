function plot_fakedatatest_CD

file_path = 'CD/faketest/';
modelname = {'1_F','1_P','1_U','2_F','2_P','2_U','3_X','3_F','3_P','3_U','4_X','4_F','4_P','4_U'};
sbjname = modelname;

nSets = 5; % number of fake data sets per model
tmp_mat = zeros(nSets,length(modelname));
LLH_map = zeros(length(modelname),length(modelname));

for i = 1:length(sbjname)
    for j = 1:nSets
        for k = 1:length(modelname)            
            fname = [file_path 'LLH_' modelname{k} '_' sbjname{i} '_' num2str(j) '.mat'];
            load(fname,'L');
            tmp_mat(j,k) = L;            
        end
    end
    tmp_mat = mean(tmp_mat);
    LLH_map(:,i) = tmp_mat - max(tmp_mat);
end

save all
qdas

% LLH_map(LLH_map<-3) = -3;
LLH_map(LLH_map<-20) = -20;

% surface plot for LLH
figure;
imagesc(LLH_map);title('model log likelihood of synthetic data (CD)');xlabel('Synthetic dataset');ylabel('Models');
