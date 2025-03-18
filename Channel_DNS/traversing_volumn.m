% traversing all combination of probes and searching for the best
% probes are placed along the central line on the end face of the volumn
% Woii junwei.chen@uc3m.es 240212
% using new datasets

%% settings
FieldFolderTrain = '../volumn_edge_new/Fields/';
FieldFolderTest  = '../volumn_edge_new/Fields_testing/';
ProbeFolderTrain = '../volumn_edge_new/Probes_linear/';
ProbeFolderTest  = '../volumn_edge_new/Probes_linear_testing/';
SaveName         = 'ProbeTraversing_volumn.mat';

S_train = 1:7200;         % snapshots for train (7200)
S_test = 1:3600;          % snapshots for test (3600)
EpL0 = 152;               % length of sensor time-series at full scale
EpB0 = 20;                % starting point of time-series at full scale
ROI = [1 88 1 88 16 27];   % first/last y, x, z (full space 88x88x44)
% x always ends at 88, z always centred at 22
NSensors = 3;             % number of sensors
ResSensors = 3;           % resolution of sensors
pTH1 = 0.7;
pTH2 = 0.7;
% row ratio at pTH1 preserves at least entries at pTH2 ratio
RModes_percentage = 95;   % energy reserved in QR-pivoting

EpL = round(EpL0/88*(ROI(4)-ROI(3)+1));
EpB = EpB0;              % probe data only exist at x = 88
tic

%% loading data
disp('loading data...');
load('geo_volumn.mat');
[X_full, Y_full, Z_full] = meshgrid(xb, yb, zb);
[X,Y,Z] = meshgrid(xb(ROI(3):ROI(4)),yb(ROI(1):ROI(2)),zb(ROI(5):ROI(6)));
u_ROI = @(u) u(ROI(1):ROI(2), ROI(3):ROI(4), ROI(5):ROI(6));

disp('loading field data for training...');
U_train = zeros(numel(X), length(S_train));
V_train = zeros(numel(X), length(S_train));
W_train = zeros(numel(X), length(S_train));
iCount = 0;
for iSnap = S_train
    iCount = iCount + 1;
    tmp = [FieldFolderTrain,'Field_',num2str(iSnap,'%.6u'),'.mat'];
    load(tmp, 'u','v','w');
    u = u_ROI(u);              U_train(:,iCount) = u(:);
    v = u_ROI(v);              V_train(:,iCount) = v(:);
    w = u_ROI(w);              W_train(:,iCount) = w(:);
end
disp('loading field data for testing...');
U_test  = zeros(numel(X), length(S_test));
V_test  = zeros(numel(X), length(S_test));
W_test  = zeros(numel(X), length(S_test));
iCount = 0;
for iSnap = S_test
    iCount = iCount + 1;
    tmp = [FieldFolderTest,'Field_',num2str(iSnap,'%.6u'),'.mat'];
    load(tmp, 'u','v','w');
    u = u_ROI(u);              U_test(:,iCount)  = u(:);
    v = u_ROI(v);              V_test(:,iCount)  = v(:);
    w = u_ROI(w);              W_test(:,iCount)  = w(:);
end
% row dir: probes, col dir: frames
disp('loading probe data for training...');
Upr_train = zeros((ROI(2)-ROI(1)+1)*EpL, length(S_train));
iCount = 0;
for iSnap = S_train
    iCount = iCount + 1;
    load([ProbeFolderTrain,'Probe_',num2str(iSnap,'%.6u'),'.mat'], 'UP');
    UP1 = UP((1:EpL)+EpB, ROI(1):ROI(2))';
    Upr_train(:,iCount) = UP1(:);
end
disp('loading probe data for testing...');
Upr_test  = zeros((ROI(2)-ROI(1)+1)*EpL, length(S_test));
iCount = 0;
for iSnap = S_test
    iCount = iCount + 1;
    load([ProbeFolderTest,'Probe_',num2str(iSnap,'%.6u'),'.mat'], 'UP');
    UP1 = UP((1:EpL)+EpB, ROI(1):ROI(2))';
    Upr_test(:,iCount)  = UP1(:);
end

clear EpB0 EpL0 FieldFol* ProbeFol* iCount iSnap tmp
clear u v w UP UP1 u_ROI
toc

%% building 2D correlation map
disp('building 2D correlation map...');
Acorr = 0.*X;
for iY = 1:size(X,1)
    for iZ = 1:size(X,3)
        tmp1 = U_train(Y==Y(iY,1,1)&Z==Z(1,1,iZ),:);
        tmp2 = U_train(X==X(1,end,1)&Y==Y(iY,1,1)&Z==Z(1,1,iZ),:);
        Acorr(iY,:,iZ) = corr(tmp1', tmp2');
    end
end
Acorr = mean(Acorr, 3);

Acorr_r = sort(Acorr, 2);
tmp = sort(Acorr_r(:, round((1-pTH2)*size(X,2))));
CorrTh = tmp(round((1-pTH1)*size(X,1)));
CorrTh = (CorrTh < 0.05)*(0.05 - CorrTh) + CorrTh;

clear Acorr_r iY iZ tmp1 tmp2 tmp
toc

%% SVD of velocity field
disp('singular value decomposition...');
global Um Vm Wm PsiU SigmaU PhiU
Um = mean(U_train, 2);   Vm = mean(V_train, 2);   Wm = mean(W_train, 2);
[PsiU, SigmaU, PhiU] = svd([U_train-Um;V_train-Vm;W_train-Wm]', 'econ');
FieldStd = std([U_train-Um;V_train-Vm;W_train-Wm], 0, 'all');
% std of the flow field
toc

%% Traversing over sensor position
disp('traversing over sensor positions...');
Row_list = 1:ResSensors:size(X,2);
AR = nchoosek(Row_list, NSensors)'; % all combinations of sensor positions
nR = size(AR, 2);
% array of error of sensor combinations of row sensors, row sensors with
% penalty from correlation value
Aerr_Row  = zeros(1, nR);
Aerr_CPen = zeros(1, nR);
Aerr_Hist = zeros(1, nR);
err_gen = @(Urecon) sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
MapShift = (round(size(X,3)/2)-1)*numel(X(:,:,1));
Map0_Row  = reshape(1:numel(X(:,:,1)), size(X(:,:,1)));
Map0_CPen = Map0_Row;
Map0_CPen(Acorr < CorrTh) = -1;
Map0_Hist = reshape(1:(ROI(2)-ROI(1)+1)*EpL, [ROI(2)-ROI(1)+1, EpL]);
warning('off');
for iProb = 1:nR
    % row sensors
    map = Map0_Row(AR(:,iProb),:) + MapShift;
    Urecon = EPOD(U_train(map(:),:), U_test(map(:),:));
    err = err_gen(Urecon);
    Aerr_Row(iProb) = err;
    fprintf([repmat('%d ',1,NSensors),'%.4f'],AR(:,iProb),err/FieldStd);
    
    % row sensors with correlation penalty
    map = Map0_CPen(AR(:,iProb),:);
    map(map == -1) = [];
    map = map + MapShift;
    Urecon = EPOD(U_train(map(:),:), U_test(map(:),:));
    err = err_gen(Urecon);
    Aerr_CPen(iProb) = err;
    fprintf(' %.4f',err/FieldStd);
    
    % sensors with time-series
    map = Map0_Hist(AR(:,iProb),:);
    Urecon = EPOD(Upr_train(map(:),:), Upr_test(map(:),:));
    err = err_gen(Urecon);
    Aerr_Hist(iProb) = err;
    fprintf(' %.4f\n',err/FieldStd);
end
warning('on');
clear iProb map Urecon err
toc

%% introducing equidistant sensors
disp('calculating result of equal distance sensors...');
tmpD = (size(X, 2) - 1)/NSensors;
RowEqi = round((tmpD/2):tmpD:size(X,2));% row selected

warning('off');
% row sensors
map = Map0_Row(RowEqi,:) + MapShift;
Urecon = EPOD(U_train(map(:),:), U_test(map(:),:));
err = err_gen(Urecon);
Derr_Row = err;
fprintf([repmat('%d ',1,NSensors),'%.4f'],RowEqi,err/FieldStd);

% row sensors with correlation penalty
map = Map0_CPen(RowEqi,:);
map(map == -1) = [];
map = map + MapShift;
Urecon = EPOD(U_train(map(:),:), U_test(map(:),:));
err = err_gen(Urecon);
Derr_CPen = err;
fprintf(' %.4f',err/FieldStd);

% sensors with time-series
map = Map0_Hist(RowEqi,:);
Urecon = EPOD(Upr_train(map(:),:), Upr_test(map(:),:));
err = err_gen(Urecon);
Derr_Hist = err;
fprintf(' %.4f\n',err/FieldStd);
warning('on');

clear tmpD map Urecon err
toc

%% error from block QR-pivoting
disp('running QR-pivoting over entire rows...');
% generating 2D slice from 3D field
iZ = 22;
sub_U_train = U_train(Z==zb(iZ),:);
sub_V_train = V_train(Z==zb(iZ),:);
sub_W_train = W_train(Z==zb(iZ),:);
sub_Um = mean(sub_U_train, 2);
sub_Vm = mean(sub_V_train, 2);
sub_Wm = mean(sub_W_train, 2);
[~, sub_SigmaU, sub_PhiU] = svd([sub_U_train-sub_Um;sub_V_train-sub_Vm;...
    sub_W_train-sub_Wm]', 'econ');
% rank in QR-pivoting
if RModes_percentage == 100
    RModes = length(sub_SigmaU);
else
    SigmaU_E = diag(sub_SigmaU.^2)/sum(sub_SigmaU.^2, 'all');
    SigmaU_EA = cumsum(SigmaU_E);
    RModes = find(SigmaU_EA <= RModes_percentage/100, 1, 'last');
    clear SigmaU_E SigmaU_EA
end
% block QR-pivoting type: sensors in the whole row,
% only for streamwise component
Phi_r = sub_PhiU(1:end/3,1:RModes);
Phi_r = RCPermutation(Phi_r, X(:,:,1));
RowQR = RowPivotQR(Phi_r*Phi_r', NSensors, size(X,1));
toc
RowQR = sort(RowQR);

warning('off');
% row sensors
map = Map0_Row(RowQR,:) + MapShift;
Urecon = EPOD(U_train(map(:),:), U_test(map(:),:));
err = err_gen(Urecon);
Perr_Row = err;
fprintf([repmat('%d ',1,NSensors),'%.4f'],RowQR,err/FieldStd);

% row sensors with correlation penalty
map = Map0_CPen(RowQR,:);
map(map == -1) = [];
map = map + MapShift;
Urecon = EPOD(U_train(map(:),:), U_test(map(:),:));
err = err_gen(Urecon);
Perr_CPen = err;
fprintf(' %.4f',err/FieldStd);

% sensors with time-series
map = Map0_Hist(RowQR,:);
Urecon = EPOD(Upr_train(map(:),:), Upr_test(map(:),:));
err = err_gen(Urecon);
Perr_Hist = err;
fprintf(' %.4f\n',err/FieldStd);
warning('on');
clear iZ sub_U_train sub_Um sub_V_train sub_Vm sub_W_train sub_Wm tmp
clear sub_SigmaU sub_PhiU Phi_r
clear map Urecon err
toc

%% save
save(SaveName, 'FieldStd', 'Acorr', 'CorrTh', 'Row_list',...
    'AR',  'Aerr_Row', 'Aerr_CPen', 'Aerr_Hist',...
    'RowEqi', 'Derr_Row', 'Derr_CPen', 'Derr_Hist',...
    'RowQR', 'Perr_Row', 'Perr_CPen', 'Perr_Hist');

%% functions

function Urecon = EPOD(ProbeDataTrain, ProbeDataTest)
% EPOD to reconstruct flow field
global Um Vm Wm PsiU SigmaU PhiU
Pm = mean(ProbeDataTrain, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeDataTrain-Pm)','econ');
Xi=PsiP'*PsiU;
Psi_recon = (ProbeDataTest-Pm)'*PhiP/SigmaP*Xi;
% [b,a] = butter(6,0.1);
% Psi_recon = filtfilt(b, a, Psi_recon);
Urecon=(Psi_recon*SigmaU*PhiU')'+[Um;Vm;Wm];
end

function A = RCPermutation(A, X)
% row and column permutation inside every column of A
% for a matrix A, every column stands for a snapshot
% where the data is in 2-D meshgrid
ColNo = (1:size(X,1):numel(X))' + ((1:size(X,1))-1);
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
