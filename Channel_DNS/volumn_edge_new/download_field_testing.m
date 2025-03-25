% 2201115
% This codes downloads the flow field for the channel flow from JHTDB
% in a volumn

HeadOfDownload;

tFolder = 'Fields_testing'; % target folder of output
% if exist(tFolder) == 7
%     rmdir(tFolder,'s');
% end
mkdir(tFolder)

dx=0.0114;
% base grid
xb=0:dx:1;
yb=0:dx:1;
zb=0:dx:0.5;
[Y_grid,X_grid,Z_grid]=meshgrid(yb,xb,zb);

% time grid
T=0:0.0065:3.25;

% space grid
Xst=0;
Zst=0;
[XstArray,ZstArray]=meshgrid(Xst,Zst);
YstArray = 0;
NSample = numel(XstArray)*length(T);

for iSample = 1:NSample
    tic
    fprintf('%d of %d samples\n',iSample,NSample);
    Xst = XstArray(mod(iSample-1,numel(XstArray))+1);
    Yst = YstArray;
    Zst = ZstArray(mod(iSample-1,numel(XstArray))+1);
    SPoints = [Xst+X_grid(:)';Yst+Y_grid(:)';Zst+Z_grid(:)'];
    t = T(floor((iSample-1)/numel(XstArray))+1);
    % t = T(ceil(rand*length(T)));
    flag = true;
    while flag
        try
            result4 =  getVelocityAndPressure (AuthKey, DataSet, t,...
                Lag6, PCHIPInt, numel(X_grid), SPoints);
            flag = false;
        end
    end
    u = reshape(result4(1,:), size(X_grid));
    v = reshape(result4(2,:), size(X_grid));
    w = reshape(result4(3,:), size(X_grid));
    p = reshape(result4(4,:), size(X_grid));
    pOutput=sprintf([tFolder,'/Field_%06d.mat'],iSample);
    save(pOutput,'u','v','w','p','xb','yb','zb',...
        'Xst','Yst','Zst','t'); 
    toc
end