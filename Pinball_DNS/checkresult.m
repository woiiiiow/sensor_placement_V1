% Woii junwei.chen@uc3m.es 221203
% check result of sensor placement

load('ProbeTraversing_R3_3pr.mat');

%% show correlation map
figure; imagesc(Acorr);
colorbar; axis equal; colormap(jet); caxis([0 1]);
hold on; contour(Acorr, (0.1:0.1:0.9), '-w', 'LineWidth', 1.5);
contour(Acorr, CorrTh:10, '-r', 'LineWidth', 1.5);
set(gca,'YDir','normal');

%% show combined scatter plot
figure;
subplot(2,2,2); histogram(Aerr_Hist/FieldStd);
xlabel('Std error from sensors with time-series');
ylabel('counts');
subplot(2,2,3); histogram(Aerr_Row/FieldStd);
xlabel('Std error from row sensors');
ylabel('counts');
subplot(2,2,1);
scatter(Aerr_Row/FieldStd,Aerr_Hist/FieldStd, 6, 'filled');
xlabel('Std error from row sensors');
ylabel('Std error from sensors with time-series');
% subplot(2,2,4);
% scatter(Aerr_Hist/FieldStd, Aerr_CPen/FieldStd, 6, 'filled');
% xlabel('Std error from sensors with time-series');
% ylabel('Std error from row sensors with correlation penalty');
subplot(2,2,4);
scatter(Aerr_CPen/FieldStd, Aerr_Hist/FieldStd, 6, 'filled');
ylabel('Std error from sensors with time-series');
xlabel('Std error from row sensors with correlation penalty');

%% show the error of equidistant sensors
subplot(2,2,1); hold on;
p1 = scatter(Derr_Row/FieldStd,Derr_Hist/FieldStd, 24, 'or','filled');
legend(p1, {'equal distant'}, 'location', 'southeast');
subplot(2,2,4); hold on;
% p1 = scatter(Derr_Hist/FieldStd,Derr_CPen/FieldStd, 24, 'or','filled');
p1 = scatter(Derr_CPen/FieldStd,Derr_Hist/FieldStd, 24, 'or','filled');
legend(p1, {'equal distant'}, 'location', 'southeast');

%% show the error of QR-pivoting
subplot(2,2,1); hold on;
p2 = scatter(Perr_Row/FieldStd,Perr_Hist/FieldStd, 24, 'og','filled');
legend([p1 p2], {'equal distant', 'QR-pivoting'}, 'location', 'southeast');
subplot(2,2,4); hold on;
% p2 = scatter(Perr_Hist/FieldStd,Perr_CPen/FieldStd, 24, 'og','filled');
p2 = scatter(Perr_CPen/FieldStd,Perr_Hist/FieldStd, 24, 'og','filled');
legend([p1 p2], {'equal distant', 'QR-pivoting'}, 'location', 'southeast');

%% show one-dimensional heat map of optimal sensor placement
NSCSe = 24; % number of sensor combination selected
table = exp(0:-1/(NSCSe-1):-1);
ResSensors = Row_list(2) - Row_list(1);

[~, iR] = mink(Aerr_Row, NSCSe);
Heat_Row = zeros(1,length(Row_list));
for i = 1:length(Row_list)
    Heat_Row(i) = sum((AR(:,iR) == Row_list(i))*table');
end
Heat_Row = interp1(Row_list, Heat_Row, 1:size(X,1), 'linear');

[~, iR] = mink(Aerr_CPen, NSCSe);
Heat_CPen = zeros(1,length(Row_list));
for i = 1:length(Row_list)
    Heat_CPen(i) = sum((AR(:,iR) == Row_list(i))*table');
end
Heat_CPen = interp1(Row_list, Heat_CPen, 1:size(X,1), 'linear');

[~, iR] = mink(Aerr_Hist, NSCSe);
Heat_Hist = zeros(1,length(Row_list));
for i = 1:length(Row_list)
    Heat_Hist(i) = sum((AR(:,iR) == Row_list(i))*table');
end
Heat_Hist = interp1(Row_list, Heat_Hist, 1:size(X,1), 'linear');

figure; hold on;
DISTMAX = 12;
DISTAMP = 1;
p1 = plot(DISTMAX - Heat_Row*DISTAMP,  1:size(X,1), 'LineWidth', 2);
p2 = plot(DISTMAX - Heat_CPen*DISTAMP, 1:size(X,1), 'LineWidth', 2);
p3 = plot(DISTMAX - Heat_Hist*DISTAMP, 1:size(X,1), 'LineWidth', 2);
legend([p1 p2 p3], {'row sensors',...
    'row sensors with correlation penalty',...
    'point sensors with time-series',...
    }, 'location', 'southwest');
set(gca, 'xTick', []);
ylim([0 77]);
xlim([0 DISTMAX]);
set(gca, 'yTick', [1:25:76]);

%% show heat map with box chart of streamwise velocity
Train_path = '../DataInterp100k/TrainData.mat';
Train = load(Train_path);
Umean = mean(Train.u,2);       Ufluc = Train.u - Umean;

% % u'/U on the x-t plane
% tmp = Ufluc./repmat(Umean,1,size(Ufluc,2));
% tmp = tmp(reshape(reshape(1:numel(X),size(X))',numel(X),1), :); % sequence
% plane_pool = reshape(tmp',size(tmp,2)*size(X,2),size(X,1));
% figure;
% boxchart(abs(plane_pool),'MarkerStyle','none','Orientation','horizontal');
% set(gca, 'xTick', []);
% set(gca, 'yTick', []);
% ylabel('box chart of |u''/U| fow every row of the field');

% u'/U by the end of every row
last_col_pool = Ufluc(end-size(X,1)+1:end,:)./...
    repmat(Umean(end-size(X,1)+1:end),1,size(Ufluc,2));
figure; hold on;
boxchart(abs(last_col_pool)',...
    'MarkerStyle','none','Orientation','horizontal');
set(gca, 'xTick', []);
set(gca, 'yTick', []);
ylabel('box chart of |u''/U| in last column');

DISTMAX = 12;
DISTAMP = 1;
p1 = plot(DISTMAX - Heat_Row*DISTAMP,  1:size(X,1), 'LineWidth', 2);
p2 = plot(DISTMAX - Heat_CPen*DISTAMP, 1:size(X,1), 'LineWidth', 2);
p3 = plot(DISTMAX - Heat_Hist*DISTAMP, 1:size(X,1), 'LineWidth', 2);
legend([p1 p2 p3], {'row sensors',...
    'row sensors with correlation penalty',...
    'point sensors with time-series',...
    }, 'location', 'southwest');

yyaxis right;
ylim([1 76]);
xlim([0 DISTMAX]);
set(gca, 'yTick', [1:15:76]);
ylabel('optimal sensers distribution');
set(gca, 'YColor', 'k')