function DataCell = load_de_delay(SbjName)

FilePath = 'de_delay/output/';
Files = dir([FilePath SbjName '_*.mat']);

% Data structure (collected by Ronald)
% data.***
% N: set size
% delay: delay time vec (all identical)
% stimvec: simulus value
% targetidx: target index
% targetval: target value
% startpos: ???
% respstartangle: orientation angle of the response probe in the beginning
% respangle: response angle value

% setting.***
% delaytime

dmat_1 = [];
dmat_2 = [];
dmat_3 = [];
dmat_6 = [];

for i = 1:length(Files)
    FileNameVec = Files(i).name;
    load([FilePath FileNameVec]);
    ErrorVec = data.respangle - data.targetval;

    % remap
    ErrorVec(ErrorVec > 90 & ErrorVec <= 180) = ErrorVec(ErrorVec > 90 & ErrorVec <= 180) - 180;
    ErrorVec(ErrorVec > -180 & ErrorVec <= -90) = ErrorVec(ErrorVec > -180 & ErrorVec <= -90) + 180;
    
    Nvec = data.N;
    DelayVec = data.delay;
    
    DelayIdx = FileNameVec(end-7);
    if DelayIdx == '1'
        dmat_1 = [dmat_1; Nvec' DelayVec' ErrorVec'];
    elseif DelayIdx == '2'
        dmat_2 = [dmat_2; Nvec' DelayVec' ErrorVec'];
    elseif DelayIdx == '3'
        dmat_3 = [dmat_3; Nvec' DelayVec' ErrorVec'];
    else
        dmat_6 = [dmat_6; Nvec' DelayVec' ErrorVec'];
    end
end

DataCellTemp = cell(1,6);

DataCellTemp{1} = dmat_1;
DataCellTemp{2} = dmat_2;
DataCellTemp{3} = dmat_3;
DataCellTemp{6} = dmat_6;

% N rows, TT columns
DataCell = cell(6,6);
% sort them by set sizes
SetSizeVec = unique(Nvec);

for T = 1:4
    if T <= 3
        TT = T;
    else
        TT = T+2;
    end
    datamat = DataCellTemp{TT};

    for n = 1:length(SetSizeVec)
        N = SetSizeVec(n);
    
        % find indices of trials with set size N
        DataCell{N,TT} = datamat(datamat(:,1)==N,:);
    end
    
end
