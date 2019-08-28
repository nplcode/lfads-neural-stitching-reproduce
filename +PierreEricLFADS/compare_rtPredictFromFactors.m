%% Load everything

PierreEricLFADS.Paper.loadPostTraining;

%%

ra = rc.findRuns('all');
D = ra.nDatasets;

thresholds = 0.1:0.05:0.95;
delta = 0.1;

%% Predict RTs from logistic regression, cross validated
clear rtPred_single;

prog = ProgressBar(D, 'Predicting RT from Factors via CIS thresholding');
for iR = 1:D
    prog.update(iR);
    r = rc.runs(iR);

    rtPred_single(iR) = PierreEricLFADS.RTPredict.run_predictRTFromFactors_CIS(r, 'thresholds', thresholds, 'delta', delta);
end
prog.finish();

rho_single = [rtPred_single.rho]';

%% Build concatenated stitched datasets

ra = rc.findRuns('all', 1);
D = ra.nDatasets;
tdSet = ra.trialDataSet;

for iD = 1:D
    tdSet{iD} = tdSet{iD}.setAxisValueList('targetDirectionName', ...
        {'Right', 'UpRight', 'Up', 'UpLeft', 'Left', 'DownLeft', 'Down', 'DownRight'});
end

%%

D = numel(tdSet);
psetsMove = cellvec(D);
psetsGo = cellvec(D);

chList = tdSet{1}.listAnalogChannelsInGroup('factors');

prog = ProgressBar(D, 'Building psets');
for iD = 1:D
    prog.update(iD);
    
    td = tdSet{iD};
    psetsGo{iD} = PopulationTrajectorySetBuilder.fromAnalogChannelsInTrialData(td, chList);
    psetsGo{iD}.minFractionTrialsForTrialAveraging = 0.8;
    
%     psetsGo{iD}.dataMean;
    
    td = tdSet{iD}.zero('Move'); % this is new
    psetsMove{iD} = PopulationTrajectorySetBuilder.fromAnalogChannelsInTrialData(td, chList);
    psetsMove{iD}.minFractionTrialsForTrialAveraging = 0.8;
    
    psetsMove{iD}.dataMean;
end
prog.finish();

%% Concatenate psets along dataset axius

psetCatMove = PopulationTrajectorySetCrossConditionUtilities.concatenateAlongNewConditionAxis(...
    psetsMove, 'dataset', {ra.datasets.name}', ...
    'equalizeTimeVectors', true, ...
    'aggregateMarks', false, 'aggregateIntervals', true);
            
% color by target using hue map
psetCatMove = psetCatMove.setConditionAppearanceFn(@LFADS_PierreExport.Appearance.colorByTarget);

% rename bases to f#
psetCatMove = psetCatMove.setBasisNamesUsingFormatString('f%d');

%% Find CIS for all move-aligned datasets

[dpca_proj, dpcaCat] = ProjDPCA_NonOrthogonal.createFromAndProject(psetCatMove, ...
    'nBasesProjPerMarginalization', [1 0 0 0], 'computeStatistics', false, 'lambda', 1e-7);

%% Concatentate go aligned

psetCatGo = PopulationTrajectorySetCrossConditionUtilities.concatenateAlongNewConditionAxis(...
    psetsGo, 'dataset', {ra.datasets.name}', ...
    'equalizeTimeVectors', true, ...
    'aggregateMarks', false, 'aggregateIntervals', true);
            
% color by target using hue map
psetCatGo = psetCatGo.setConditionAppearanceFn(@LFADS_PierreExport.Appearance.colorByTarget);

% rename bases to f#
psetCatGo = psetCatGo.setBasisNamesUsingFormatString('f%d');

dpcaCatGo = dpca_proj.projectPopulationTrajectorySet(psetCatGo);

%% Do per dataset rt predicition

nThresh = numel(thresholds);

rho_stitched_mat = nan(D, nThresh);

prog = ProgressBar(D, 'Predicting RT from Factors via CIS thresholding');
for iR = 1:D
    prog.update(iR);
    
    td = rc.runs(iR).trialDataSet{1};
    rt = td.getParam('rt');

    dpca = dpca_proj.projectPopulationTrajectorySet(psetsGo{iR}, 'computeStatistics', false);
        
    data = dpca.dataByTrial{1};
    tvec = dpca.tvecDataByTrial{1};

    % try to flip it 
    meanCIS = mean(data, 1);
    if mean(meanCIS(tvec < 0)) > mean(meanCIS(tvec > 0))
        data = -data;
        meanCIS = -meanCIS;
    end

    tMask = tvec > -50 & tvec < 400;
    if iR == 1
        lo = quantile(meanCIS(tMask), 0.1);
        hi = quantile(meanCIS(tMask), 0.9);
    end
    data = (data - lo) ./ (hi - lo);
    
    rtPredCell = cellvec(numel(thresholds));
    rho_vals = nanvec(numel(thresholds));

    for ithr = 1:numel(thresholds)
        rtPredCell{ithr} = TrialDataUtilities.Data.findThresholdCrossingsLowThenHigh(data', tvec, thresholds(ithr), thresholds(ithr)+delta)';
        mask = ~isnan(rt) & ~isnan(rtPredCell{ithr});
        if mean(~isnan(rtPredCell{ithr})) < 0.7
            continue;
        end

        rho_stitched_mat(iR, ithr) = corr(rtPredCell{ithr}(mask), rt(mask));
    end
   
%     [rho, idx] = nanmax(rho_vals);
%     rtPred = rtPredCell{idx};
%     
%     clf;
%     pt(2, tvec, data, 'coloreval', rt);
%     hold on;
%     horzLine(thresholds(idx));
%     hold off;
% 
%     rtPred = TrialDataUtilities.Data.findThresholdCrossingsLowThenHigh(data', tvec, 0.6, 0.8)';
%     
%     rho_stitched(iR) = nancorr(rtPred, rt);
end
prog.finish();

[~, idxThresh] = max(nansum(rho_stitched_mat, 1));
rho_stitched = rho_stitched_mat(:, idxThresh);
thresh_stitched = thresholds(idxThresh);

%% Save rho values

out_file = fullfile(results_dir, 'rtPredict_rho_from_CIS.mat');
save(out_file, 'rho_stitched', 'thresh_stitched', 'rho_single');

%% Load results
out_file = fullfile(results_dir, 'rtPredict_rho_from_CIS.mat');
load(out_file);


%% Plot rho comparison

figure(1); clf
scatter(rho_single, rho_stitched, 15, 'filled', 'MarkerEdgeColor', 'w', 'MarkerEdgeAlpha', 0.5, 'LineWidth', 0.1);
set(gca, 'TickDir', 'out');
% title(sprintf('Comparison ', opts.lag, opts.binWidth, opts.source), 'Background', 'none');
title('RT Predictions');
xlabel('Single session \rho')
ylabel('Stitched \rho');
xlim([0 1]);
ylim([0 1]);
hold on;
hl = TrialDataUtilities.Plotting.identityLine;
uistack(hl, 'bottom');
% niceGrid;
axis equal;
xlim([0 1]);
ylim([0 1]);
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
ax = AutoAxis.replace;
figSize([5 5]);
hold off;
grid on

%% Save fig

saveFigure(fullfile(fig_dir, 'rtPrediction_fromCIS_rho_compare_single_vs_stitching'), 'painters', true, 'upsample', 4);

%% Print result

fname = fullfile(results_dir, 'rt_results.txt');
fid = fopen(fname, 'w');
fprintf(fid, 'Mean improvement: %.4g\n', mean(rho_stitched - rho_single));
fprintf(fid, 'p value signrank: %.8g\n', signrank(rho_single, rho_stitched));
fclose(fid);

dbtype(fname);

%% Plot single dataset example

iR = 12; % fixed
% [rho_best, iR] = max(rho_stitched);
rho_this = rho_stitched(iR);

ds = dc.datasets(12);
dsName = sprintf('%s%s', ds.subject(1), datestr(ds.datenum, 'YYYYmmdd'));

td = ra.trialDataSet{iR};
rt = td.getParam('rt');

dpca = dpca_proj.projectPopulationTrajectorySet(psetsGo{iR}, 'computeStatistics', false);

data = dpca.dataByTrial{1};
tvec = dpca.tvecDataByTrial{1};

% try to flip it 
meanCIS = mean(data, 1);
if mean(meanCIS(tvec < 0)) > mean(meanCIS(tvec > 0))
    data = -data;
    meanCIS = -meanCIS;
end
tMask = tvec > -50 & tvec < 400;
data = (data - lo) ./ (hi - lo);
rtPred = TrialDataUtilities.Data.findThresholdCrossingsLowThenHigh(data', tvec, thresh_stitched, thresh_stitched+0.1)';

cmap = flipud(TrialDataUtilities.Color.cmocean('haline', 300));
clim = [220 520];

figure(2); clf;
PierreEricLFADS.Utils.pt(2, tvec, data, 'colormap', cmap, 'colorlim', clim, 'coloreval', rt, 'alpha', 0.5);

title(sprintf('Dataset %s CIS, \\rho = %0.2f', dsName, rho_this));

% hc = colorbar;
% ylabel(hc, 'RT (ms)');
% hc.Box = 'off';

hold on
horzLine(thresh_stitched, 'Color', [0.3 0.3 0.3]);
hold off;
xlim([-500 700]);
ylim([-0.2 1.2]);

%%
xlabel('Time from Go Cue (ms)');
ylabel('CIS (a.u.)');
set(gca, 'YTick', []);
aa = AutoAxis.replace();

aa.addColorbar('cmap', cmap, 'orientation', 'vertical', 'width', 0.2, 'units', 'ms', 'labelCenter', 'RT');
aa.axisMarginRight = 1.5;
grid on;
aa.update();
aa.update(); % bug in auto axis code

figSize([7 5]);

%% Save fig

saveFigure(fullfile(fig_dir, 'rtPrediction_exampleDataset'), 'painters', true);


