% Woii junwei.chen@uc3m.es 221203
% check result of sensor placement

% load('ProbeTraversing_R3_3pr_edge_NM.mat');
load('ProbeTraversing_volumn_Z44.mat');

%% show correlation map
figure; imagesc(Acorr);
colorbar; axis equal; colormap(jet); caxis([0 1]);
hold on; contour(Acorr, (0.1:0.1:0.9), '-w', 'LineWidth', 1.5);
contour(Acorr, CorrTh:10, '-r', 'LineWidth', 1.5);
% set(gca,'YDir','normal');

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
subplot(2,2,4);
scatter(Aerr_Hist/FieldStd, Aerr_CPen/FieldStd, 6, 'filled');
xlabel('Std error from sensors with time-series');
ylabel('Std error from row sensors with correlation penalty');

%% show the error of equidistant sensors
subplot(2,2,1); hold on;
p1 = scatter(Derr_Row/FieldStd,Derr_Hist/FieldStd, 24, 'or','filled');
legend(p1, {'equal distant'}, 'location', 'southeast');
subplot(2,2,4); hold on;
p1 = scatter(Derr_Hist/FieldStd,Derr_CPen/FieldStd, 24, 'or','filled');
legend(p1, {'equal distant'}, 'location', 'southeast');

%% show the error of QR-pivoting
subplot(2,2,1); hold on;
p2 = scatter(Perr_Row/FieldStd,Perr_Hist/FieldStd, 24, 'og','filled');
legend([p1 p2], {'equal distant', 'QR-pivoting'}, 'location', 'southeast');
subplot(2,2,4); hold on;
p2 = scatter(Perr_Hist/FieldStd,Perr_CPen/FieldStd, 24, 'og','filled');
legend([p1 p2], {'equal distant', 'QR-pivoting'}, 'location', 'southeast');

%% show one-dimensional heat map of optimal sensor placement
NSCSe = 24; % number of sensor combination selected
Height = size(Acorr,1);
table = exp(0:-1/(NSCSe-1):-1);
ResSensors = Row_list(2) - Row_list(1);

[~, iR] = mink(Aerr_Row, NSCSe);
Heat_Row = zeros(1,length(Row_list));
for i = 1:length(Row_list)
    Heat_Row(i) = sum((AR(:,iR) == Row_list(i))*table');
end
Heat_Row = interp1(Row_list, Heat_Row, 1:Height, 'linear');

[~, iR] = mink(Aerr_CPen, NSCSe);
Heat_CPen = zeros(1,length(Row_list));
for i = 1:length(Row_list)
    Heat_CPen(i) = sum((AR(:,iR) == Row_list(i))*table');
end
Heat_CPen = interp1(Row_list, Heat_CPen, 1:Height, 'linear');

[~, iR] = mink(Aerr_Hist, NSCSe);
Heat_Hist = zeros(1,length(Row_list));
for i = 1:length(Row_list)
    Heat_Hist(i) = sum((AR(:,iR) == Row_list(i))*table');
end
Heat_Hist = interp1(Row_list, Heat_Hist, 1:Height, 'linear');

figure; hold on;
DISTMAX = 0.3;
DISTAMP = 1/50;
p1 = plot(DISTMAX - Heat_Row*DISTAMP,  1:Height, 'LineWidth', 2);
p2 = plot(DISTMAX - Heat_CPen*DISTAMP, 1:Height, 'LineWidth', 2);
p3 = plot(DISTMAX - Heat_Hist*DISTAMP, 1:Height, 'LineWidth', 2);
legend([p1 p2 p3], {'row sensors',...
    'row sensors with correlation penalty',...
    'point sensors with time-series',...
    }, 'location', 'northwest');
set(gca, 'xTick', []);
ylim([0 89]);
xlim([0 DISTMAX]);
set(gca, 'yTick', [1:29:88]);
set(gca,'YDir','reverse');

