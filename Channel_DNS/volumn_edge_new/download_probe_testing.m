                                                                                                                                           % 2201107
% This codes downloads the virtual probe series for the channel flow
% from JHTDB for a volumn

tFolder = 'Probes_linear_testing'; % target folder of output
mkdir(tFolder);

% flow field data
fTemplate = 'Fields_testing/F*';
Af = dir(fTemplate);

HeadOfDownload;

% frame number of each probe records
EpisodeLength = 200;
EpisodeRetro = 20;
dt = 0.0065;

% base grid
dx=0.0114;
x_grid = 1;
y_grid = 0:dx:1;
z_grid = 0.2394;
[X_pr,Y_pr,Z_pr]=meshgrid(x_grid,y_grid,z_grid);

for iErgo = 1:7200
    tic
    load([Af(iErgo).folder,'/',Af(iErgo).name], 't','Xst','Yst','Zst');
    fprintf('%d of %d samples\n',iErgo,length(Af));
    SPoints = [Xst+X_pr(:)';Yst+Y_pr(:)';Zst+Z_pr(:)'];
    UP = zeros(EpisodeLength, numel(X_pr)); VP = UP; WP = VP;
    tp = t - EpisodeRetro*dt;
    fprintf('downloading probe data at the frame behind the snapshot: ');
    for iEp = 1:EpisodeLength
        fprintf('%03d', iEp-1-EpisodeRetro);
        flag = true;
        while flag
            try
                tmp_velo = getVelocity...
                    (AuthKey,DataSet,tp,Lag6,PCHIPInt,numel(X_pr),SPoints);
                flag = false;
            end
        end
        UP(iEp,:) = reshape(tmp_velo(1,:), [1, numel(X_pr)]);
        VP(iEp,:) = reshape(tmp_velo(2,:), [1, numel(X_pr)]);
        WP(iEp,:) = reshape(tmp_velo(3,:), [1, numel(X_pr)]);
        tp = tp + dt;
        fprintf('\b\b\b');
    end
    fprintf('\n');
    
    pOutput = [tFolder,'/',strrep(Af(iErgo).name,'Field','Probe')];
    save(pOutput,'dt','x_grid','y_grid','z_grid','Xst','Yst','Zst','t',...
        'EpisodeLength', 'EpisodeRetro','UP','VP','WP');
    toc
end
