% woii 240109 junwei.chen@uc3m.es

% settings
FieldFolder = '../database/';
SubFieldFolder = {'PIV21/', 'PIV22/', 'PIV23/', 'PIV24/'};
FileName = 'flow.mat';
N_train = 4800;                % size of train set
TestFieldFolder = 'PIV24/';    % folder of test set
S_test = 1001:1500;            % test set
ROI = [24 108 24 73 2 8];     % subdomain of the PIV field
% [first-last row, first-last column first-last slice]
X_pr = ROI(4);
Y_pr = ROI(1):ROI(2);
Z_pr = 5;
EpL = round(80*(ROI(4)-ROI(3))/108);
% EpisodeLength, length of probe segment
NSensors = 3; 
ResSensors = 3;
pTH1 = 0.7;
pTH2 = 0.7; % row ratio at pTH1 preserves at least entries at pTH2 ratio
RModes_percentage = 95;        % energy reserved in QR-pivoting
SaveName = sprintf('ProbeTraversing_R%d_%dpr_NM.mat', ResSensors,NSensors);
tic

% generating list of training data
sub_nf = zeros(1, length(SubFieldFolder));
cand_train = [];
for iF = 1:length(SubFieldFolder)
    if iF == 1
        load([FieldFolder, SubFieldFolder{iF}, FileName], 'X','Y','Z',...
            'Mm_per_px_ratio','Sample_rate','Vector_Spacing');
    end
    tmp = load([FieldFolder, SubFieldFolder{iF}, FileName], 'AFrame');
    sub_nf(iF) = length(tmp.AFrame);
    sub_list = 1:sub_nf(iF);
    sub_list(end-EpL+1:end) = 0;
    if strcmp(SubFieldFolder{iF}, TestFieldFolder)
        sub_list(S_test) = 0;
    end
    cand_train = [cand_train, sub_list(sub_list~=0)+sum(sub_nf(1:iF-1))];
end
NFrame = sum(sub_nf);
S_train = cand_train(sort(randperm(length(cand_train), N_train)));

%% loading data sets
disp('loading data...');
map_field = false(size(X));
map_field(ROI(1):ROI(2), ROI(3):ROI(4), ROI(5):ROI(6)) = true;
map_probe = false(size(X));
map_probe(Y_pr, X_pr, Z_pr) = true;
U_train = zeros(sum(map_field, 'all'), N_train);
V_train = U_train;             W_train = U_train;
Upr_train = zeros(sum(map_probe, 'all')*EpL, N_train);
Upr_test  = zeros(sum(map_probe, 'all')*EpL, length(S_test));
iflag     = 0;
for iF = 1:length(SubFieldFolder)
    fprintf('loading block %d...\n', iF);
    tmp = load([FieldFolder, SubFieldFolder{iF}, FileName], 'U','V','W');
    itable = S_train > sum(sub_nf(1:iF-1)) & S_train < sum(sub_nf(1:iF))+1;
    sub_index = S_train(itable)-sum(sub_nf(1:iF-1));
    U_train(:, iflag+(1:sum(itable))) = tmp.U(map_field, sub_index);
    V_train(:, iflag+(1:sum(itable))) = tmp.V(map_field, sub_index);
    W_train(:, iflag+(1:sum(itable))) = tmp.W(map_field, sub_index);
    for isub = 1:sum(itable)
        isnap = S_train(iflag+isub) - sum(sub_nf(1:iF-1));
        tmp2 = tmp.U(map_probe, isnap-1+(1:EpL));
        Upr_train(:,iflag+isub) = tmp2(:);
    end
    iflag     = iflag     + sum(itable);
    if strcmp(SubFieldFolder{iF}, TestFieldFolder)
        U_test = tmp.U(map_field, S_test);
        V_test = tmp.V(map_field, S_test);
        W_test = tmp.W(map_field, S_test);
        for isub = 1:length(S_test)
            tmp2 = tmp.U(map_probe, S_test(isub)-1+(1:EpL));
            Upr_test(:,isub) = tmp2(:);
        end
    end
end
clear cand_train FieldFolder FileName i* tmp* sub_index sub_list
toc

%% building 2D correlation map
disp('building 2D correlation map...');
Acorr = zeros(ROI(2)-ROI(1)+1,ROI(4)-ROI(3)+1);
X_sub = X(ROI(1):ROI(2), ROI(3):ROI(4), ROI(5):ROI(6));
Y_sub = Y(ROI(1):ROI(2), ROI(3):ROI(4), ROI(5):ROI(6));
Z_sub = Z(ROI(1):ROI(2), ROI(3):ROI(4), ROI(5):ROI(6));
X = X_sub;   Y = Y_sub;   Z = Z_sub;

for iY = 1:size(Acorr,1)
    iZ = 5;
    tmp1 = U_train(Y==Y(iY,1,1)&Z==Z(1,1,iZ),:);
    tmp2 = U_train(X==X(1,end,1)&Y==Y(iY,1,1)&Z==Z(1,1,iZ),:);
    Acorr(iY,:) = corr(tmp1', tmp2');
end

% Corrth: correlation value under this value is to be mapped out
Acorr_r = sort(Acorr, 2);
tmp = sort(sort(Acorr_r(:, round((1-pTH2)*size(X,2)))));
CorrTh = tmp(round((1-pTH1)*size(X,1)));
CorrTh = (CorrTh < 0.05)*(0.05 - CorrTh) + CorrTh;

figure; hold on; imagesc(Acorr);
axis equal; colorbar; colormap(jet); caxis([0 1]);
contour(Acorr, CorrTh:10, '-r', 'LineWidth', 1.5);

clear iY iZ tmp1 tmp2 tmp
toc

%% SVD of velocity field
disp('singular value decomposition...');
Um = mean(U_train, 2);
Vm = mean(V_train, 2);
Wm = mean(W_train, 2);
[PsiU, SigmaU, PhiU] = svd([U_train-Um;V_train-Vm;W_train-Wm]', 'econ');
% std of the flow field
FieldStd = std([U_train-Um;V_train-Vm;W_train-Wm], 0, 'all');
toc

%% Traversing over sensor position
disp('traversing over sensor positions...');
Row_list = 1:ResSensors:size(X,1);
AR = nchoosek(Row_list, NSensors)'; % all combinations of sensor positions
nR = size(AR, 2);
% array of error of sensor combinations of row sensors, row sensors with
% penalty from correlation value
Aerr_Row  = zeros(1, nR);
Aerr_CPen = zeros(1, nR);
Aerr_Hist = zeros(1, nR);
MapShift  = (round(size(X,3)/2)-1)*numel(X(:,:,1));
Map0_Row  = reshape(1:numel(X(:,:,1)), size(X(:,:,1)));
Map0_CPen = Map0_Row;
Map0_CPen(Acorr < CorrTh) = -1;
Map0_Hist = reshape(1:size(Y,1)*EpL, [size(Y,1), EpL]);
warning('off');
for iProb = 1:nR
    % row sensors
    % sensor positions in the domain
    map = Map0_Row(AR(:,iProb),:) + MapShift;
    % train of EPOD
    ProbeData = U_train(map(:),:);
    Pm = mean(ProbeData, 2);
    [PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
    Xi=PsiP'*PsiU;
    % test of EPOD
    ProbeData = U_test(map(:),:);
    Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm;Wm];
    err = sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
    Aerr_Row(iProb) = err;
    fprintf([repmat('%d ',1,NSensors),'%.4f'],AR(:,iProb),err/FieldStd);
    
    % sensors with correlation penalty
    % sensor positions in the domain
    map = Map0_CPen(AR(:,iProb),:);
    map(map == -1) = [];
    map = map + MapShift;
    % train of EPOD
    ProbeData = U_train(map(:),:);
    Pm = mean(ProbeData, 2);
    [PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
    Xi=PsiP'*PsiU;
    % test of EPOD
    ProbeData = U_test(map(:),:);
    Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm;Wm];
    err = sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
    Aerr_CPen(iProb) = err;
    fprintf(' %.4f',err/FieldStd);
    
    % sensors with time-series
    map = Map0_Hist(AR(:,iProb),:);
    ProbeData = Upr_train(map(:),:);
    Pm = mean(ProbeData, 2);
    [PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
    Xi=PsiP'*PsiU;
    % test of EPOD
    ProbeData = Upr_test(map(:),:);
    Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm;Wm];
    err = sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
    Aerr_Hist(iProb) = err;
    fprintf(' %.4f\n',err/FieldStd);
end
warning('on');
clear iProb map Pm PsiP SigmaP PhiP ProbeData Urecon err
toc

%% introducing equidistant sensors
disp('calculating result of equal distance sensors...');
tmpD = (size(X, 1) - 1)/NSensors;
RowEqi = round((tmpD/2):tmpD:size(X,1));% row selected

warning('off');
% row sensors
map = Map0_Row(RowEqi,:) + MapShift; % sensor positions in the domain
% train of EPOD
ProbeData = U_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = U_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm;Wm];
err = sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
Derr_Row = err;
fprintf([repmat('%d ',1,NSensors),'%.4f'],RowEqi,err/FieldStd);
% sensors with correlation penalty
map = Map0_CPen(RowEqi,:); % sensor positions in the domain
map(map == -1) = [];
map = map + MapShift;
% train of EPOD
ProbeData = U_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = U_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm;Wm];
err = sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
Derr_CPen = err;
fprintf(' %.4f',err/FieldStd);
% sensors with time-series
map = Map0_Hist(RowEqi,:);
ProbeData = Upr_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = Upr_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm;Wm];
err = sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
Derr_Hist = err;
fprintf(' %.4f\n',err/FieldStd);
warning('on');
clear tmpD map Pm PsiP SigmaP PhiP ProbeData Urecon err
toc

%% error from block QR-pivoting
disp('running QR-pivoting over entire rows...');
% generating 2D slice from 3D field
iZ = 5;
sub_U_train = U_train(Z==Z(1,1,iZ),:);
sub_V_train = V_train(Z==Z(1,1,iZ),:);
sub_W_train = W_train(Z==Z(1,1,iZ),:);
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
RowQR = RowPivotQR(Phi_r*Phi_r', NSensors, size(X,2));
toc
RowQR = sort(RowQR);

warning('off');
% row sensors
map = Map0_Row(RowQR,:) + MapShift; % sensor positions in the domain
% train of EPOD
ProbeData = U_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = U_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm;Wm];
err = sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
Perr_Row = err;
fprintf([repmat('%d ',1,NSensors),'%.4f'],RowQR,err/FieldStd);
% sensors with correlation penalty
map = Map0_CPen(RowQR,:); % sensor positions in the domain
map(map == -1) = [];
map = map + MapShift;
% train of EPOD
ProbeData = U_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = U_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm;Wm];
err = sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
Perr_CPen = err;
fprintf(' %.4f',err/FieldStd);
% sensors with time-series
map = Map0_Hist(RowQR,:);
ProbeData = Upr_train(map(:),:);
Pm = mean(ProbeData, 2);
[PsiP, SigmaP, PhiP] = svd((ProbeData-Pm)','econ');
Xi=PsiP'*PsiU;
% test of EPOD
ProbeData = Upr_test(map(:),:);
Urecon=((ProbeData-Pm)'*PhiP/SigmaP*Xi*SigmaU*PhiU')'+[Um;Vm;Wm];
err = sqrt(mean(([U_test;V_test;W_test]-Urecon).^2, 'all'));
Perr_Hist = err;
fprintf(' %.4f\n',err/FieldStd);
warning('on');
clear iZ sub_U_train sub_Um sub_V_train sub_Vm sub_W_train sub_Wm tmp
clear sub_SigmaU sub_PhiU Phi_r
clear map Pm PsiP SigmaP PhiP ProbeData Urecon err tmp
toc

%% save
save(SaveName, 'FieldStd', 'Acorr', 'CorrTh', 'Row_list', 'ROI',...
    'AR',  'Aerr_Row', 'Aerr_CPen', 'Aerr_Hist',...
    'RowEqi', 'Derr_Row', 'Derr_CPen', 'Derr_Hist',...
    'RowQR', 'Perr_Row', 'Perr_CPen', 'Perr_Hist');

%% functions

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

