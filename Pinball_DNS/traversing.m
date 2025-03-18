% traversing all combination of probes and searching for the best
% Woii junwei.chen@uc3m.es 221202

%% settings
Train_path = 'TrainData.mat';
Test_path  = 'TestData.mat';
Probe_path = 'ProbeHistory_everyrow.mat';

S_train = -1;           % snapshots for train, -1 means all
S_test = -1;            % snapshots for test, -1 means all
NSensors = 3;           % number of sensors
ResSensors = 3;         % resolution of sensor position
EpisodeLength = 80;     % length of sensor time-series
RModes_percentage = 95; % energy reserved in QR-pivoting
pTH1 = 0.7;
pTH2 = 0.7; % row ratio at pTH1 preserves at least entries at pTH2 ratio
tic

%% loading data
disp('loading data...');
Train = load(Train_path);
Test  = load(Test_path);
Probe = load(Probe_path);
X = Train.X;                   Y = Train.Y;
Test.Nimg(end-100:end)= [];
Test.u(:,end-100:end) = [];
Test.v(:,end-100:end) = [];
Test.p(:,end-100:end) = [];
U_train = Train.u;             V_train = Train.v;
U_test = Test.u;               V_test = Test.v;
Upr_train = zeros(EpisodeLength*size(X,1), length(Train.Nimg));
Vpr_train = Upr_train;
for iC = 1:length(Train.Nimg)
    iP = find(Probe.Nimg == Train.Nimg(iC),1)-1;
    tmp = Probe.u(:,iP + (1:EpisodeLength));
    Upr_train(:,iC) = tmp(:);
end
Upr_test = zeros(EpisodeLength*size(X,1), length(Test.Nimg));
Vpr_test = Upr_test;
for iC = 1:length(Test.Nimg)
    iP = find(Probe.Nimg == Test.Nimg(iC),1)-1;
    tmp = Probe.u(:,iP + (1:EpisodeLength));
    Upr_test(:,iC) = tmp(:);
end

if S_train(1) ~= -1
    U_train = U_train(:,S_train);
    V_train = V_train(:,S_train);
    Upr_train = Upr_train(:,S_train);
end
if S_test(1) ~= -1
    U_test = U_test(:,S_test);
    V_test = V_test(:,S_test);
    Upr_test = Upr_test(:,S_test);
end

toc
clear *_path iC iP tmp
clear Train Test Probe S_train S_test

%% building correlation map
disp('building correlation map...');

% only stream wise velocity
Acorr = 0.*X;
IndexTable = 1:size(X,1):numel(X); % the index of first row of X in U_train
for iRow = 1:size(X,1)
    Acorr(iRow,:) = corr(U_train(IndexTable+iRow-1,:)',...
        U_train(numel(X)-size(X,1)+iRow,:)');
end
toc

% % average over streamwise and crosswise velocity
% Acorr1 = 0.*X;
% Acorr2 = 0.*X;
% IndexTable = 1:size(X,1):numel(X); % the index of first row of X in U_train
% for iRow = 1:size(X,1)
%     Acorr1(iRow,:) = corr(U_train(IndexTable+iRow-1,:)',...
%         U_train(numel(X)-size(X,1)+iRow,:)');
%     Acorr2(iRow,:) = corr(V_train(IndexTable+iRow-1,:)',...
%         V_train(numel(X)-size(X,1)+iRow,:)');
% end
% Acorr = (Acorr1 + Acorr2)./2;
% toc

% figure; imagesc(Acorr);
% colorbar; axis equal; colormap(jet); caxis([0 1]);
% hold on; contour(Acorr, (0.1:0.1:0.9), '-w');

%% regular of correlation penalty
% Corrth: correlation value under this value is to be mapped out
Acorr_r = sort(Acorr, 2);
tmp = sort(sort(Acorr_r(:, round((1-pTH2)*size(X,2)))));
CorrTh = tmp(round((1-pTH1)*size(X,1)));
CorrTh = (CorrTh < 0.05)*(0.05 - CorrTh) + CorrTh;

%% SVD of velocity field
disp('singular value decomposition...');
Um = mean(U_train, 2);         Vm = mean(V_train, 2);
[PsiU, SigmaU, PhiU] = svd([U_train-Um;V_train-Vm]', 'econ');
FieldStd = std([U_train-Um;V_train-Vm], 0, 'all'); % std of the flow field
toc

%% Traversing over sensor position
disp('traversing over sensor positions...');
Row_list = 1:ResSensors:size(X,1);
AR = nchoosek(Row_list, NSensors)'; % all combinations of sensor positions
nR = size(AR, 2);
% array of error of sensor combinations of row sensors, row sensors with
% penalty from correlation value and point sensors with time-series
Aerr_Row  = zeros(1, nR);
Aerr_CPen = zeros(1, nR);
Aerr_Hist = zeros(1, nR);
Map0_Row  = reshape(1:numel(X), size(X));
Map0_CPen = Map0_Row;
Map0_Hist = reshape(1:size(X,1)*EpisodeLength, [size(X,1), EpisodeLength]);
Map0_CPen(Acorr < CorrTh) = -1;
warning('off');
for iProb = 1:nR
    % row sensors
    map = Map0_Row(AR(:,iProb),:); % sensor positions in the domain
    % train of EPOD
    ProbeData = U_train(map(:),:);
    Pm = mean(ProbeData, 2);
    [PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
    Xi=PsiP'*PsiU;
    % test of EPOD
    ProbeData = U_test(map(:),:);
    Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm];
    err = sqrt(mean(([U_test;V_test]-Urecon).^2, 'all'));
    Aerr_Row(iProb) = err;
    fprintf([repmat('%d ',1,NSensors),'%.4f'],AR(:,iProb),err/FieldStd);
    
    % sensors with correlation penalty
    map = Map0_CPen(AR(:,iProb),:);
    map(map == -1) = [];
    % train of EPOD
    ProbeData = U_train(map(:),:);
    Pm = mean(ProbeData, 2);
    [PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
    Xi=PsiP'*PsiU;
    % test of EPOD
    ProbeData = U_test(map(:),:);
    Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm];
    err = sqrt(mean(([U_test;V_test]-Urecon).^2, 'all'));
    Aerr_CPen(iProb) = err;
    fprintf(' %.4f',err/FieldStd);
    
    % sensors with time-series
    map = Map0_Hist(AR(:,iProb),:);
    % train of EPOD
    ProbeData = Upr_train(map(:),:);
    Pm = mean(ProbeData, 2);
    [PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
    Xi=PsiP'*PsiU;
    % test of EPOD
    ProbeData = Upr_test(map(:),:);
    Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm];
    err = sqrt(mean(([U_test;V_test]-Urecon).^2, 'all'));
    Aerr_Hist(iProb) = err;
    fprintf(' %.4f\n',err/FieldStd);

    % toc
end
warning('on');
toc

%% introducing equidistant sensors
disp('calculating result of equal distance sensors...');
tmpD = (size(X, 1) - 1)/NSensors;
RowEqi = round((tmpD/2):tmpD:size(X,1));% row selected

warning('off');
% row sensors
map = Map0_Row(RowEqi,:); % sensor positions in the domain
% train of EPOD
ProbeData = U_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = U_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm];
err = sqrt(mean(([U_test;V_test]-Urecon).^2, 'all'));
Derr_Row = err;
fprintf([repmat('%d ',1,NSensors),'%.4f'],RowEqi,err/FieldStd);
% sensors with correlation penalty
map = Map0_CPen(RowEqi,:);
map(map == -1) = [];
% train of EPOD
ProbeData = U_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = U_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm];
err = sqrt(mean(([U_test;V_test]-Urecon).^2, 'all'));
Derr_CPen = err;
fprintf(' %.4f',err/FieldStd);
% sensors with time-series
map = Map0_Hist(RowEqi,:);
% train of EPOD
ProbeData = Upr_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = Upr_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm];
err = sqrt(mean(([U_test;V_test]-Urecon).^2, 'all'));
Derr_Hist = err;
fprintf(' %.4f\n',err/FieldStd);
warning('on');
toc


%% error from QR-pivoting
disp('running QR-pivoting over entire rows...');
% rank in QR-pivoting
if RModes_percentage == 100
    RModes = length(SigmaU);
else
    SigmaU_E = diag(SigmaU.^2)/sum(SigmaU.^2, 'all');
    SigmaU_EA = cumsum(SigmaU_E);
    RModes = find(SigmaU_EA <= RModes_percentage/100, 1, 'last');
    clear SigmaU_E SigmaU_EA
end
% sensors in the whole row, only for streamwise component
Phi_r = PhiU(1:end/2,1:RModes);
Phi_r = RCPermutation(Phi_r, X);
RowQR = RowPivotQR(Phi_r*Phi_r', NSensors, size(X,2));
toc
RowQR = sort(RowQR);
warning('off');
% row sensors
map = Map0_Row(RowQR,:); % sensor positions in the domain
% train of EPOD
ProbeData = U_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = U_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm];
err = sqrt(mean(([U_test;V_test]-Urecon).^2, 'all'));
Perr_Row = err;
fprintf([repmat('%d ',1,NSensors),'%.4f'],RowQR,err/FieldStd);
% sensors with correlation penalty
map = Map0_CPen(RowQR,:);
map(map == -1) = [];
% train of EPOD
ProbeData = U_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = U_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm];
err = sqrt(mean(([U_test;V_test]-Urecon).^2, 'all'));
Perr_CPen = err;
fprintf(' %.4f',err/FieldStd);
% sensors with time-series
map = Map0_Hist(RowQR,:);
% train of EPOD
ProbeData = Upr_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = Upr_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm];
err = sqrt(mean(([U_test;V_test]-Urecon).^2, 'all'));
Perr_Hist = err;
fprintf(' %.4f\n',err/FieldStd);
warning('on');
toc

%% saving result
disp('saving result...')
save(sprintf('ProbeTraversing_R%d_%dpr_NM.mat', ResSensors, NSensors),...
    'X', 'Y', 'Acorr', 'CorrTh', 'FieldStd', ...
    'Row_list', 'AR', 'Aerr_Row', 'Aerr_CPen', 'Aerr_Hist',...
    'RowEqi', 'Derr_Row', 'Derr_CPen', 'Derr_Hist',...
    'RowQR', 'Perr_Row', 'Perr_CPen', 'Perr_Hist');
toc

%% functions

function A = RCPermutation(A, X)
% row and column permutation inside every column of A
% for a matrix A, every column stands for a snapshot
% where the data is in 2-D meshgrid
ColNo = [1:size(X,1):numel(X)]' + ([1:size(X,1)]-1);
A = A(ColNo(:),:);
end

function plist = RowPivotQR(A, NSensors, LenRow)
% self version of QR Pivoting, every row of original data set is regarded
% as a group
% input: A - m x n matrix, m <= n; NSensors - number of sensors;
% LenRow - row length of the field.
% output: plist 1 x NSensors vector, row selected in QR Pivoting
pivot = 1:size(A,2)/LenRow; % sequence of rows in pivoting
plist = 1:NSensors;
for k = 1:NSensors
    % calculating the mean-square of blocks in A, each block stands for a
    % row in QR pivoting, but since A is normalized during pivoting,
    % calcalation only starts from kth block in bottom-right corner.
    tmp = sum(A((k-1)*LenRow+1:end,(k-1)*LenRow+1:end).^2);
    tmp = sum(reshape(tmp, [LenRow,numel(tmp)/LenRow]));
    % finding the block of maximum mean-square value
    p = find(tmp==max(tmp),1)-1+k;
    % swapping the kth and pth block in A
    tmp = A(:,(k-1)*LenRow+[1:LenRow]);
    A(:,(k-1)*LenRow+[1:LenRow]) = A(:,(p-1)*LenRow+[1:LenRow]);
    A(:,(p-1)*LenRow+[1:LenRow]) = tmp;
    % swapping the sequence in pivot
    tmp = pivot(p);   pivot(p) = pivot(k);   pivot(k) = tmp;
    % recording the sequence of block with maximum mean-square in current
    % loop into plist
    plist(k) = pivot(k);
    % Householder reflections for kth block in A after swapping, so that
    % the lower triangular part of the block are transformed to be 0
    for ksub = (k-1)*LenRow+[1:LenRow]
        u = A(ksub:end,ksub);
        u = u - [norm(u);zeros(size(A,1)-ksub,1)];
        u = u./norm(u);
        Asub = u'*A(ksub:end,ksub:end);
        A(ksub:end,ksub:end) = A(ksub:end,ksub:end) - 2*u*Asub;
    end
end
end
