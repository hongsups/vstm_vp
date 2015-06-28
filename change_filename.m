clear all

% filepath = 'CD/faketest/';
% 
% model_group_name = {'1' '2' '3' '4'};
% for i = 1:4
%     files = dir([filepath 'LLH_' model_group_name{i} '*.mat']);
%     for j = 1:length(files)
%         nameOld = files(j).name;
%         nameNew = strrep(nameOld,['LLH_' model_group_name{i}],'LLH');
%         movefile([filepath nameOld], [filepath nameNew]);
%     end
% end

% model_group_name = {'IL' 'SA' 'EP' 'VP'};
% for i = 1:4
%     files = dir([filepath 'LLH_' model_group_name{i} '*.mat']);
%     for j = 1:length(files)
%         nameOld = files(j).name;
%         if i == 1
%             nameNew = strrep(nameOld,['LLH_' model_group_name{i}],'LLH_1');
%             nameNewNew = strrep(nameNew,'_3S','_');
%             movefile([filepath nameOld], [filepath nameNewNew]);
%         else
%             nameNew = strrep(nameOld,['LLH_' model_group_name{i}],'LLH');
%             movefile([filepath nameOld], [filepath nameNew]);
%         end
%     end
% end

% filepath = 'conjunction/table/';
% 
% iVec = 1:15;
% jVec = 1:15;
% table_N = 'T_VPconj_r1r2_';
% 
% for i = 1:length(iVec)
%     for j = 1:length(jVec)
%         files = dir([filepath table_N num2str(i) '_' num2str(j) '_' '*.mat']);
%         nameOld = files.name;
%         nameNew = strrep(nameOld,'_500.mat','.mat');
%         movefile([filepath nameOld], [filepath nameNew]);
%     end
% end

% filepath = 'conjunction/table/';
% 

PathName = 'de_delay/';
DelayVec = [1 2 3 6];
Nvec = [1 2 4 6];
SbjNameCell = {'HJK','HS','JYP','RRS','YS'};

for i = 1:5
    SbjName = SbjNameCell{i};
    for j = 1:4
    DelayIdx = DelayVec(j);
        for k = 1:4
            NIdx = Nvec(k);
            FileName_VP = [PathName 'de_delay_mix_' SbjName '_' num2str(DelayIdx*1000) '_' num2str(NIdx)];
            fileVP1 = dir([FileName_VP '.mat']);
            nameOld1 = fileVP1.name;
            nameNew1 = strrep(nameOld1, 'mix','mix_np');
            movefile([PathName nameOld1],[PathName nameNew1]);
            
            fileVP2 = dir([FileName_VP '.fig']);
            nameOld2 = fileVP2.name;
            nameNew2 = strrep(nameOld2, 'mix','mix_np');
            movefile([PathName nameOld2],[PathName nameNew2]);
        end
    end
end

