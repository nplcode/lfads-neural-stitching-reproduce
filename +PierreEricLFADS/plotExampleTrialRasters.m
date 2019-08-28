PierreEricLFADS.Paper.loadPostTraining;

ra = rc.findRuns('all');
D = ra.nDatasets;
tdSet = ra.trialDataSet;

%% Plot all rasters on same figure

figure(1); clf;
figSize([15 9]);
R = 6;
C = 7;
g = AutoAxisGrid(gcf, R, C);
M = R * C;
for iR = 1:M
    ax = g.axisAt(iR);
    td = tdSet{iR};
    td = td.filter('targetDirectionName', 'Up');
    td = td.mark('GoCue', 'as', 'Go').mark('Move');
    td.plotSingleTrialRaster('validTrialIdx', 1, 'quick', true);
    
    ax.Units = 'centimeters';
    ax.LooseInset = [0.2 0.2 0.2 0.2];
    ax.XTick = [];
    ax.YTick = [];
    delete(ax.XLabel);
    ax.XRuler.Visible = 'off';
    ax.YRuler.Visible = 'off';
    
    if iR >= M-C+1
        td.setupTimeAxis('showScaleBar', iR == D, 'showRanges', false);
        aa = AutoAxis(ax);
        aa.axisMargin = [0.1 1 0.1 0.1];
    end
end

aa = AutoAxis(ax);
aa.addAutoScaleBarX('length', 200, 'units', 'ms');
aa.axisMargin = [0.1 1 0.1 0.1];
AutoAxis.updateFigure();

drawnow;
g.updatePositions();

%%

fig_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper_v5';
saveFigure(fullfile(fig_dir, 'exampleRasters_UpReach'), 'painters', true);

%% Plot 5 random rasters on same figure

figure(1); clf;
figSize([14 3]);
N = 5;
g = AutoAxisGrid(gcf, 1, N);

rng(1);
idx = randsample(45, N);

for i = 1:N
    iR = idx(i);
    
    ax = g.axisAt(i);
    td = tdSet{iR};
%     td = td.filter('targetDirectionName', 'Up');
    rt = td.getParam('rt');
    td = td.withTrials(rt > 300).filter('targetDirectionName', 'Up');
    td = td.mark('GoCue').mark('Move');
    td.plotSingleTrialRaster('validTrialIdx', 1, 'quick', true);
    
    ax.Units = 'centimeters';
    ax.LooseInset = [0.2 0.2 0.2 0.2];
    ax.XTick = [];
    ax.YTick = [];
    delete(ax.XLabel);
    ax.XRuler.Visible = 'off';
    ax.YRuler.Visible = 'off';

    td.selectValidTrials(1).setupTimeAxis('showScaleBar', iR == D, 'showRanges', false);
    aa = AutoAxis(ax);
    aa.axisMargin = [0.1 1 0.1 0.1];
end

aa = AutoAxis(ax);
aa.addAutoScaleBarX('length', 200, 'units', 'ms');
aa.axisMargin = [0.1 1 0.3 0.1];
AutoAxis.updateFigure();

drawnow;
g.updatePositions();

%%

fig_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper';
saveFigure(fullfile(fig_dir, 'exampleRastersSubset_UpReach'), 'painters', true);

%% Plot vertically stacked rasters with ellipsis

figure(1); clf;
figSize([4 6]);
N = 5;
g = AutoAxisGrid(gcf, N, 1);

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

idx = [1:N-1, D];

for i = 1:N
    iR = idx(i);
    
    ax = g.axisAt(i);
    td = tdSet{iR};
%     td = td.filter('targetDirectionName', 'Up');
    rt = td.getParam('rt');
    td = td.withTrials(rt > 300).filter('targetDirectionName', 'Up');
    td = td.mark('GoCue').mark('Move');
    td.plotSingleTrialRaster('validTrialIdx', 1, 'quick', true);
    
    ax.Units = 'centimeters';
    ax.LooseInset = [0.2 0.2 0.2 0.2];
    ax.XTick = [];
    ax.YTick = [];
    delete(ax.XLabel);
    ax.XRuler.Visible = 'off';
    ax.YRuler.Visible = 'off';
    
    ylabel(sprintf('Sess. %d', iR));
    
    aa = AutoAxis(ax);
    if iR == D
        % last one
        ht = title('...');
        td.selectValidTrials(1).setupTimeAxis('showScaleBar', iR == D, 'showRanges', false, 'showLabels', true, 'style', 'marker');
    else
        td.selectValidTrials(1).setupTimeAxis('showScaleBar', iR == D, 'showRanges', false, 'showLabels', false, 'style', 'marker');
        aa.axisMargin = [0.3 0.2 0.1 0.05];
    end
end

aa = AutoAxis(ax);
aa.addAutoScaleBarX('length', 200, 'units', 'ms');
aa.axisMargin = [0.2 0.5 0.2 0.3];
aa.axisPaddingBottom = 0.1;
AutoAxis.updateFigure();

drawnow;
g.updatePositions();

ht.Units = 'centimeters';
ht.Position(2) = ht.Position(2) + 0.05;
ht.FontWeight = 'bold';

%%
saveFigure(fullfile(fig_dir, 'exampleRastersVertical_UpReach'), 'painters', true);


%% Plot factors

figure(1); clf;
figSize([15 15]);
g = AutoAxisGrid(gcf, 9, 5);
prog = ProgressBar(D, 'Plotting sessions');
for iR = 1:D
    prog.update(iR);
    ax = g.axisAt(iR);
    td = tdSet{iR};
    td = td.mark('GoCue').mark('Move');
    rt = td.getParam('rt');
    td = td.withTrials(rt > 300).filter('targetDirectionName', 'Up');
%     td.plotSingleTrialAnalogChannelGroup('factors', 'validTrialIdx', 1, ...
%         'commonTime', true, 'timeDelta', 20, 'quick', true, ...
%         'intercalate', true, 'normalize', false, 'colormapStacked', @TrialDataUtilities.Color.hslmap);
    
    % T x N
    [mat, tvec] = td.selectValidTrials(1).getAnalogChannelGroupAsTensor('factors', 'timeDelta', 20);

    mat = squeeze(mat);
    pmat(mat', 'x', tvec, 'addColorbar', false);
    
    ax.Units = 'centimeters';
    ax.LooseInset = [0.2 0.2 0.2 0.2];
    ax.XTick = [];
    ax.YTick = [];
    delete(ax.XLabel);
    ax.XRuler.Visible = 'off';
    ax.YRuler.Visible = 'off';
    
    if iR >= D-4
        td.setupTimeAxis('showScaleBar', iR == D);
        aa = AutoAxis(ax);
        aa.axisMargin = [0.1 1 0.1 0.1];
    end

end
prog.finish();

% equalize color limits
clim = cell2mat(cellfun(@(ax) caxis(ax), g.handles(:), 'UniformOutput', false));
clim = [min(clim(:, 1)) max(clim(:, 2))];
cellfun(@(ax) caxis(ax, clim), g.handles);

aa = AutoAxis(ax);
aa.addAutoScaleBarX('length', 200, 'units', 'ms');
aa.axisMargin = [0.1 1 0.1 0.1];

AutoAxis.updateFigure();

drawnow;
g.updatePositions();

%%

fig_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper';
saveFigure(fullfile(fig_dir, 'exampleFactors_UpReach'), 'painters', true);

%% Plot Rates

%% Plot factors

figure(1); clf;
figSize([15 15]);
g = AutoAxisGrid(gcf, 9, 5);
prog = ProgressBar(D, 'Plotting sessions');
for iR = 1:D
    prog.update(iR);
    ax = g.axisAt(iR);
    td = tdSet{iR};
    td = td.mark('GoCue').mark('Move');
    td = td.filter('targetDirectionName', 'Up');
%     td.plotSingleTrialAnalogChannelGroup('factors', 'validTrialIdx', 1, ...
%         'commonTime', true, 'timeDelta', 20, 'quick', true, ...
%         'intercalate', true, 'normalize', false, 'colormapStacked', @TrialDataUtilities.Color.hslmap);
    
    % T x N
    [mat, tvec] = td.selectValidTrials(1).getAnalogChannelGroupAsTensor('rates', 'timeDelta', 20);

    mat = squeeze(mat);
    pmat(mat', 'x', tvec, 'addColorbar', false);
    
    ax.Units = 'centimeters';
    ax.LooseInset = [0.2 0.2 0.2 0.2];
    ax.XTick = [];
    ax.YTick = [];
    delete(ax.XLabel);
    ax.XRuler.Visible = 'off';
    ax.YRuler.Visible = 'off';
    
    if iR >= D-4
        td.setupTimeAxis('showScaleBar', iR == D);
        aa = AutoAxis(ax);
        aa.axisMargin = [0.1 1 0.1 0.1];
    end

end
prog.finish();

% equalize color limits
% clim = cell2mat(cellfun(@(ax) caxis(ax), g.handles(:), 'UniformOutput', false));
% clim = [min(clim(:, 1)) max(clim(:, 2))];
% cellfun(@(ax) caxis(ax, clim), g.handles);

aa = AutoAxis(ax);
aa.addAutoScaleBarX('length', 200, 'units', 'ms');
aa.axisMargin = [0.1 1 0.1 0.1];

AutoAxis.updateFigure();

drawnow;
g.updatePositions();

%%

fig_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper';
saveFigure(fullfile(fig_dir, 'exampleRates_UpReach'), 'painters', true);