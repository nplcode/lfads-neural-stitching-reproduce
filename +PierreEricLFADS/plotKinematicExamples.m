%% load everything
PierreEricLFADS.Paper.loadPostTraining;

%%

ra = rc.findRuns('all');
D = ra.nDatasets;

%% Pick the dataset that is representative


out_file = fullfile(results_dir, 'lfads_kinematics_fromFactors.mat');
r2_lfads = load(out_file);

r2_stitched_mean = mean(r2_lfads.r2_stitched_factors, 2);
r2_single_mean = mean(r2_lfads.r2_single_factors, 2);
delta = r2_stitched_mean - r2_single_mean;

% mask = abs(delta - median(delta)) < 0.1;
mask = delta > 0.2;
score = r2_stitched_mean;
score(~mask) = -Inf;
[val, iDS] = max(score);


tbl = table(r2_single_mean, r2_stitched_mean, delta, 'RowNames', arrayfun(@num2str, 1:44, 'UniformOutput', false))

% fixing at dataset 32
% (0.58 single, 0.789 stitched, 0.20 delta)
iDS = 27;

fprintf('Choosing dataset %d with stitched = %.2f, delta = %.2f\n', iDS, r2_stitched_mean(iDS), delta(iDS));
% Selection was iDS == 31 
% Switching to iDS == 15

% iDS = 15;

out_file = fullfile(results_dir, 'gpfa_decode_r2only.mat');
data = load(out_file);

r2_gpfa_mean = mean(data.r2_gpfa, 2);

out_file = fullfile(results_dir, 'smoothed_neural_decode_r2only.mat');
data = load(out_file);

r2_smoothed_mean = mean(data.r2_smoothed, 2);

mat = cat(2, r2_smoothed_mean, r2_gpfa_mean, r2_single_mean, r2_stitched_mean);

r2_this = mat(iDS, :);

%% Load the kinematic predictions for everything

ra = rc.runs(end);
out_path = fullfile(results_dir, 'kinematic_predictions');

out_file = fullfile(out_path, sprintf('kinematics_true_%s.mat', ra.datasets(iDS).name));
tmp = load(out_file);
T_true = tmp.T_true;

out_file = fullfile(out_path, sprintf('kinematics_lfadsSingleFromFactors_%s.mat', ra.datasets(iDS).name));
tmp = load(out_file);
T_lfads_single = tmp.T_lfads_single;
save(out_file, 'T_lfads_single');

out_file = fullfile(out_path, sprintf('kinematics_lfadsStitchedFromFactors_%s.mat', ra.datasets(iDS).name));
tmp = load(out_file);
T_lfads_stitched = tmp.T_lfads_stitched;

out_file = fullfile(out_path, sprintf('kinematics_gpfa_%s.mat', ra.datasets(iDS).name));
tmp = load(out_file);
T_gpfa = tmp.T_gpfa;
    
out_file = fullfile(out_path, sprintf('kinematics_smoothed_neural_%s.mat', ra.datasets(iDS).name));
tmp = load(out_file);
T_smoothed_neural = tmp.T_smoothed_neural;

binWidth = 20;

%% Plot all conditions

td = ra.trialDataSet{iDS};
td = LFADS_PierreExport.Condition.groupCenterOut(td);
conditionIdx = td.conditionIdx;
cmap = td.conditionColors;

Tset = {T_true; T_smoothed_neural; T_gpfa; T_lfads_single; T_lfads_stitched};
% decode_names = {'Hand', 'Smoothed', 'GPFA', sprintf('Single session\nLFADS'), sprintf('Multi-session\nLFADS')};
decode_names = {'Hand', sprintf('Smoothed\nr^2 = %.2f', r2_this(1)), ...
    sprintf('GPFA\nr^2 = %.2f', r2_this(2)), ...
    sprintf('Single Session LFADS\nr^2 = %.2f', r2_this(3)), ...
    sprintf('Stitched LFADS\nr^2 = %.2f', r2_this(4))};
nT = numel(Tset);

figure(3); clf;
figSize([13.5 3]);

g = AutoAxisGrid(1, nT);

for t = 1:nT
    ax = g.axisAt(t);
    PierreEricLFADS.Utils.plot_kinematic_decode_from_tstruct(Tset{t}, binWidth, ...
        'colorByCondition', true, 'conditionIdx', conditionIdx, 'colorMap', cmap, 'alpha', 0.7);
    box off;
    xlabel(decode_names{t});
    axis tight;
    axis equal;
    ax.XTick = []; ax.XRuler.Axle.Visible = 'off';
    ax.YTick = []; ax.YRuler.Axle.Visible = 'off';
    ax.LooseInset = ax.TightInset;
    
    if t < nT
        ax = AutoAxis();
        ax.axisMargin = [0 0.7 0 0];
    else
        ax  = AutoAxis.replaceScaleBars(ax, 'mm', 'mm');
        ax.axisMargin = [0 0.7 0.5 0];
        ax.axisPaddingRight = 0.35;
        ax.axisPaddingBottom = 0;
        ax.scaleBarLenX = 40;
        ax.scaleBarLenY = 40;
    end
end

hax = cat(1, g.handles{:});
lims = [cell2mat(xlim(hax)), cell2mat(ylim(hax))];
% L = max(abs(lims(:)));
L = 120;
xlim(hax, [-L L]);
ylim(hax, [-L L]);

% ax.addXLabel('Multi-session');
AutoAxis.updateFigure();
g.updatePositions();

%%

saveFigure(fullfile(fig_dir, sprintf('exampleKinematics_%s', dc.datasets(iDS).name)), 'painters', true);


%% Plot all datasets

prog = ProgressBar(D, 'Plotting datasets');
for iDS = 1:D
    prog.update(iDS);
    
    % Load the kinematic predictions for everything

    out_file = fullfile(out_path, sprintf('kinematics_true_%s.mat', ra.datasets(iDS).name));
    tmp = load(out_file);
    T_true = tmp.T_true;

    out_file = fullfile(out_path, sprintf('kinematics_lfadsSingleFromFactors_%s.mat', ra.datasets(iDS).name));
    tmp = load(out_file);
    T_lfads_single = tmp.T_lfads_single;
    save(out_file, 'T_lfads_single');

    out_file = fullfile(out_path, sprintf('kinematics_lfadsStitchedFromFactors_%s.mat', ra.datasets(iDS).name));
    tmp = load(out_file);
    T_lfads_stitched = tmp.T_lfads_stitched;

    out_file = fullfile(out_path, sprintf('kinematics_gpfa_%s.mat', ra.datasets(iDS).name));
    tmp = load(out_file);
    T_gpfa = tmp.T_gpfa;

    out_file = fullfile(out_path, sprintf('kinematics_smoothed_neural_%s.mat', ra.datasets(iDS).name));
    tmp = load(out_file);
    T_smoothed_neural = tmp.T_smoothed_neural;

    td = ra.trialDataSet{iDS};
    td = LFADS_PierreExport.Condition.groupCenterOut(td);
    conditionIdx = td.conditionIdx;
    cmap = td.conditionColors;

    Tset = {T_true; T_smoothed_neural; T_gpfa; T_lfads_single; T_lfads_stitched};
    decode_names = {'Hand', 'Smoothed', 'GPFA', sprintf('LFADS'), sprintf('Multi LFADS')};
    nT = numel(Tset);

    figure(1); clf;
    figSize([7.5*2.54/2 10*2.54/11]);

    g = AutoAxisGrid(1, nT);

    for t = 1:nT
        ax = g.axisAt(t);
        PierreEricLFADS.Utils.plot_kinematic_decode_from_tstruct(Tset{t}, binWidth, ...
            'colorByCondition', true, 'conditionIdx', conditionIdx, 'colorMap', cmap);
        box off;
        xlabel(decode_names{t});
        axis tight;
        axis equal;
        ax.XTick = []; ax.XRuler.Axle.Visible = 'off';
        ax.YTick = []; ax.YRuler.Axle.Visible = 'off';
        ax.LooseInset = ax.TightInset;

        if t == 1
            ax = AutoAxis();
            ax.ylabel(sprintf('Session %d', iDS));
            hy = get(gca, 'YLabel');
            hy.FontWeight = 'bold';
            ax.axisMargin = [0.3 0.4 0 0];
        elseif t < nT
            ax = AutoAxis();
            ax.axisMargin = [0 0.4 0 0];
        else
            ax  = AutoAxis();
            ax.addScaleBarY('units', 'mm', 'length', 40);
            ax.axisMargin = [0 0.4 0.4 0];
            ax.scaleBarLenX = 40;
            ax.scaleBarLenY = 40;
        end
    end

    hax = cat(1, g.handles{:});
    lims = [cell2mat(xlim(hax)), cell2mat(ylim(hax))];
    % L = max(abs(lims(:)));
    L = 120;
    xlim(hax, [-L L]);
    ylim(hax, [-L L]);

    % ax.addXLabel('Multi-session');
    AutoAxis.updateFigure();
    g.updatePositions();

    fig_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper_v5/all_kinematic_examples';
    saveFigure(fullfile(fig_dir, sprintf('kin_%02d', iDS)), 'ext', {'pdf'}, 'painters', true);
end
prog.finish();