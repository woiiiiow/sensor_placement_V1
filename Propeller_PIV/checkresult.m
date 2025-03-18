% Woii junwei.chen@uc3m.es 230503
% check result of sensor placement

load('ProbeTraversing_R3_3pr_NM_2D.mat');
% load('geo_volumn.mat');

%% show correlation map

figure; imagesc(Acorr);
colorbar; axis equal; colormap(jet); caxis([0 1]);
hold on; contour(Acorr, (0.1:0.1:0.9), '-w', 'LineWidth', 1.5);
contour(Acorr, CorrTh:10, '-r', 'LineWidth', 1.5);
set(gca,'YDir','normal');

%% show scatter plot
figure;
histogram(Aerr_Row/FieldStd);
xlabel('Std error from row sensors');
ylabel('counts');
figure;
histogram(Aerr_CPen/FieldStd);
xlabel('Std error from row sensors masked');
ylabel('counts');
figure;
histogram(Aerr_Hist/FieldStd);
xlabel('Std error from sensors with time-series');
ylabel('counts');

figure; hold on;
scatter(Aerr_Row/FieldStd,Aerr_Hist/FieldStd, 6, 'filled');
xlabel('Std error from row sensors');
ylabel('Std error from sensors with time-series');
figure; hold on;
scatter(Aerr_CPen/FieldStd,Aerr_Hist/FieldStd, 6, 'filled');
xlabel('Std error from row sensors masked');
ylabel('Std error from sensors with time-series');
figure; hold on;
scatter(Aerr_Row/FieldStd,Aerr_CPen/FieldStd, 6, 'filled');
xlabel('Std error from row sensors');
ylabel('Std error from row sensors masked');

%% show scatter plot
figure; hold on;
scatter(Aerr_Row/FieldStd,Aerr_Hist/FieldStd, 6, 'filled');
xlabel('Std error from row sensors');
ylabel('Std error from probe with time-series');
p1 = scatter(Derr_Row/FieldStd,Derr_Hist/FieldStd, 24, 'or','filled');
p2 = scatter(Perr_Row/FieldStd,Perr_Hist/FieldStd, 24, 'og','filled');
legend([p1 p2], {'equal distant', 'QR-pivoting'}, 'location', 'southeast');

figure; hold on;
scatter(Aerr_CPen/FieldStd,Aerr_Hist/FieldStd, 6, 'filled');
xlabel('Std error from row sensors masked');
ylabel('Std error from probe with time-series');
p1 = scatter(Derr_CPen/FieldStd,Derr_Hist/FieldStd, 24, 'or','filled');
p2 = scatter(Perr_CPen/FieldStd,Perr_Hist/FieldStd, 24, 'og','filled');
legend([p1 p2], {'equal distant', 'QR-pivoting'}, 'location', 'southeast');

%% show one-dimensional heat map of optimal sensor placement
NSCSe = 24; % number of sensor combination selected
table = exp(0:-1/(NSCSe-1):-1);
ResSensors = Row_list(2) - Row_list(1);
HeightList = 1:ROI(2)-ROI(1)+1;

[~, iR] = mink(Aerr_Row, NSCSe);
Heat_Row = zeros(1,length(Row_list));
for i = 1:length(Row_list)
    Heat_Row(i) = sum((AR(:,iR) == Row_list(i))*table');
end
Heat_Row = interp1(Row_list, Heat_Row, HeightList, 'linear');

[~, iR] = mink(Aerr_CPen, NSCSe);
Heat_CPen = zeros(1,length(Row_list));
for i = 1:length(Row_list)
    Heat_CPen(i) = sum((AR(:,iR) == Row_list(i))*table');
end
Heat_CPen = interp1(Row_list, Heat_CPen, HeightList, 'linear');

[~, iR] = mink(Aerr_Hist, NSCSe);
Heat_Hist = zeros(1,length(Row_list));
for i = 1:length(Row_list)
    Heat_Hist(i) = sum((AR(:,iR) == Row_list(i))*table');
end
Heat_Hist = interp1(Row_list, Heat_Hist, HeightList, 'linear');

figure; hold on;
DISTMAX = 120;
DISTAMP = 12;
p1 = plot(DISTMAX - Heat_Row*DISTAMP,  HeightList, 'LineWidth', 2);
p2 = plot(DISTMAX - Heat_CPen*DISTAMP, HeightList, 'LineWidth', 2);
p3 = plot(DISTMAX - Heat_Hist*DISTAMP, HeightList, 'LineWidth', 2);
legend([p1 p2 p3], {'row sensors',...
    'row sensors with correlation penalty',...
    'point sensors with time-series',...
    }, 'location', 'southwest');
set(gca, 'xTick', []);
ylim([0 ROI(2)-ROI(1)+2]);
xlim([0 DISTMAX]);
set(gca, 'yTick', '');