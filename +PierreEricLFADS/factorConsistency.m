%% Load everything

PierreEricLFADS.Paper.loadPostTraining;

%% Go Cue full period

ra = rc.findRuns('all', 1);
D = ra.nDatasets;
tdSet = ra.trialDataSet;

% for iD = 1:D
%     tdSet{iD} = LFADS_PierreExport.Condition.groupCenterOut(tdSet{iD});
% end

%% build Move aligned psets

D = numel(tdSet);
psetsMove = cellvec(D);

chList = tdSet{1}.listAnalogChannelsInGroup('factors');

prog = ProgressBar(D, 'Building psets');
for iD = 1:D
    prog.update(iD);

    td = tdSet{iD}.zero('Move');
    psetsMove{iD} = PopulationTrajectorySetBuilder.fromAnalogChannelsInTrialData(td, chList);
    psetsMove{iD}.minFractionTrialsForTrialAveraging = 0.8;
    
    psetsMove{iD}.dataMean;
end
prog.finish();

%% Concatenate psets along dataset axis

psetCatMove = PopulationTrajectorySetCrossConditionUtilities.concatenateAlongNewConditionAxis(...
    psetsMove, 'dataset', {ra.datasets.name}', ...
    'equalizeTimeVectors', true, ...
    'aggregateMarks', false, 'aggregateIntervals', false);
            
% color by target using hue map
psetCatMove = psetCatMove.setConditionAppearanceFn(@LFADS_PierreExport.Appearance.colorByTarget);

% rename bases to f#
psetCatMove = psetCatMove.setBasisNamesUsingFormatString('f%d');

% psetCatMoveShort = psetCatMove.manualSliceOrExpandTimeWindow(-200, 300);
% [psetCatMoveZ, tr_zscore] = psetCatMove.zscoreByBasis();
% [pca_proj, pcaCatMove] = ProjPCA.createFromAndProject(psetCatMoveShort, 'computeStatistics', false);

[pca_proj, pcaCatMove] = ProjPCA.createFromAndProject(psetCatMove, 'computeStatistics', false);


%%
% pcaMoveViaGo = pca_proj_go.projectPopulationTrajectorySet(psetCatMove);
% clf;
% pcaMoveViaGo.plotBases('basisIdx', 1:3);
% figure(1);
%%
clf;
pcaCatMove.plotStateSpace('useThreeVector', false, 'alpha', 0.7, 'timeDelta', 20, 'LineWidth', 0.5);
view(-29, 9);
return;

%%
clf;
pcaCatMove.plotStateSpace('useThreeVector', false, 'alpha', 0.7, 'LineWidth', 0.5);
view(-135, -10);
tv = ThreeVector();
tv.vectorLength = 0.3;
figSize([5 5]);

%%
clf
pcaCatMove.plotStateSpace('LineWidth', 0.5, 'startShowOnData', true, 'stopShowOnData', true, 'zeroShowOnData', true, 'useThreeVector', false);              

%%

[dpca_proj, dpcaCatMove] = ProjDPCA_NonOrthogonal.createFromAndProject(ps, tCatMove, 'computeStatistics', false, 'nBasesProjPerMarginalization', [2 2 0 0]);

%% saveFigure

fig_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper_v6';
saveFigure(fullfile(fig_dir, 'factorTrajectoryPlot'), 'painters', true);

%% Go cue aligned
D = numel(tdSet);
psetsGo = cellvec(D);

chList = tdSet{1}.listAnalogChannelsInGroup('factors');

prog = ProgressBar(D, 'Building psets');
for iD = 1:D
    prog.update(iD);

    td = tdSet{iD};
    psetsGo{iD} = PopulationTrajectorySetBuilder.fromAnalogChannelsInTrialData(td, chList);
    psetsGo{iD}.dataMean;
end
prog.finish();

%% Concatenate psets along dataset axis

psetCatGo = PopulationTrajectorySetCrossConditionUtilities.concatenateAlongNewConditionAxis(...
    psetsGo, 'dataset', {ra.datasets.name}', ...
    'equalizeTimeVectors', true, ...
    'aggregateMarks', false, 'aggregateIntervals', false);
            
% color by target using hue map
psetCatGo = psetCatGo.setConditionAppearanceFn(@LFADS_PierreExport.Appearance.colorByTarget);

% rename bases to f#
psetCatGo = psetCatGo.setBasisNamesUsingFormatString('f%d');

[pca_proj_go, pcaCatGo] = ProjPCA.createFromAndProject(psetCatGo.zscoreByBasis, 'computeStatistics', false);

%%
cla
pcaCatGo.manualSliceOrExpandTimeWindow(0, 700).plotStateSpace('LineWidth', 0.5, 'markSize', 20, 'startShowOnData', true, 'stopShowOnData', true, 'zeroShowOnData', true, 'useThreeVector', false);              



%%
clf;
pcaCatGo.plotStateSpace('useThreeVector', false, 'alpha', 0.7, 'timeDelta', 20, 'LineWidth', 0.5);

%%

clf;
pcaSlicedGo = pcaCatGo.manualSliceOrExpandTimeWindow(0, 490);
pcaSlicedGo.plotStateSpace( 'alpha', 0.7, 'timeDelta', 20, 'LineWidth', 0.5);

%% export factor data for chethan, trial averaged

nDirections = 8;
nDatasets = D;
tensor = pcaCatMove.dataMean{1};
sz = size(tensor);
factors_pc_by_direction_by_dataset_by_time = reshape(tensor(1:3, :, :), [3 nDirections nDatasets sz(3)]);

info = LFADS_PierreExport.Appearance.getColorSet();
colors = info.centerOutMap;
targets = info.centerOutDirectionList;
datasets = {dc.datasets.name}';

save(fullfile(fig_dir, 'stitched_pca_factor_export.mat'), ...
    'factors_pc_by_direction_by_dataset_by_time', 'colors', 'targets', 'datasets');

%% export PCA'd single trials

D = numel(tdSet);
pcaMoves = cellvec(D);
pca_proj_trunc = pca_proj.truncateOutputBases(3);

prog = ProgressBar(D, 'Building psets');
for iD = 1:D
    prog.update(iD);
    pcaMoves{iD} = pca_proj_trunc.projectPopulationTrajectorySet(psetsMove{iD});
end
prog.finish();

%%
conditionIdx_byDataset = cellfun(@(pset) pset.conditionIdxByTrial{1}, psetsMove, 'UniformOutput', false);
info = LFADS_PierreExport.Appearance.getColorSet();
colors = info.centerOutMap;
targets = info.centerOutDirectionList;
datasets = {dc.datasets.name}';

moveAligned_byDataset_trials_by_time_by_factorPC = cellfun(@(pset) cat(3, pset.dataByTrial{:}), pcaMoves, 'UniformOutput', false);
moveAligned_byDataset_time_vector = cellfun(@(pset) pset.tvecDataByTrial, pcaMoves, 'UniformOutput', false);

save(fullfile(fig_dir, 'single_trial_factorPC_export.mat'), ...
    'moveAligned_byDataset_trials_by_time_by_factorPC', 'moveAligned_byDataset_time_vector', ...
    'conditionIdx_byDataset', 'colors', 'targets', 'datasets');

%% export spike trains and kinematics

[spikesByDataset_trialByChannel, kinematicsByDataset_trialsByTimeByChannel, kinematicsByDataset_time, moveOnsetByDataset] = cellvec(D);
kinematicsChannels = {'handX', 'handY', 'handVelocityX', 'handVelocityY'};

for iD = 1:D
    td = tdSet{iD};
    spikesByDataset_trialByChannel{iD} = td.getSpikeTimes(td.listSpikeChannels);
    [kinematicsByDataset_trialsByTimeByChannel{iD}, kinematicsByDataset_time{iD}] = ...
        td.getAnalogMultiAsTensor(kinematicsChannels, 'timeDelta', 20);
    
    moveOnsetByDataset{iD} = td.getEventFirst('Move');
end
%%

save(fullfile(fig_dir, 'moveOnsetTimes_export.mat'), 'moveOnsetByDataset');

save(fullfile(fig_dir, 'spikes_kinematics_export.mat'), ...
    'spikesByDataset_trialByChannel', 'kinematicsByDataset_trialsByTimeByChannel', ...
    'kinematicsByDataset_time', ...
    'conditionIdx_byDataset', 'colors', 'targets', 'datasets');

%% export single trials with condition labels

conditionIdx_byDataset = cellfun(@(pset) pset.conditionIdxByTrial{1}, psetsGo, 'UniformOutput', false);

gcAligned_byDataset_trials_by_time_by_factor = cellfun(@(pset) cat(3, pset.dataByTrial{:}), psetsGo, 'UniformOutput', false);
gcAligned_byDataset_time_vector = cellfun(@(pset) pset.tvecDataByTrial, psetsGo, 'UniformOutput', false);

moveAligned_byDataset_trials_by_time_by_factor = cellfun(@(pset) cat(3, pset.dataByTrial{:}), psetsMove, 'UniformOutput', false);
moveAligned_byDataset_time_vector = cellfun(@(pset) pset.tvecDataByTrial, psetsMove, 'UniformOutput', false);

pcaDecoderKbyN = pca_proj.decoderKbyN(1:3, :);

save(fullfile(fig_dir, 'single_trial_factor_export.mat'), ...
    'gcAligned_byDataset_trials_by_time_by_factor', 'gcAligned_byDataset_time_vector', ...
    'moveAligned_byDataset_trials_by_time_by_factor', 'moveAligned_byDataset_time_vector', ...
    'pcaDecoderKbyN', 'conditionIdx_byDataset', 'colors', 'targets', 'datasets');

%% export single trial


%%
[dpca_proj, dpcaCatMove] = ProjDPCA_NonOrthogonal.createFromAndProject(psetCatMove, 'lambda', 1e-7, 'computeStatistics', false, 'nBasesProjPerMarginalization', 3, 'lambda');

