

%%


fig_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper_v11_resubmit';
mkdirRecursive(fig_dir);
results_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper_v11_resubmit/results';
mkdirRecursive(results_dir);

%% Load R^2 for smoothed neural, gpfa, lfads single and stitched

out_file = fullfile(results_dir, 'lfads_kinematics_fromFactors.mat');
data = load(out_file);

r2_stitched_mean = mean(data.r2_stitched_factors, 2);
r2_single_mean = mean(data.r2_single_factors, 2);

out_file = fullfile(results_dir, 'gpfa_decode_r2only.mat');
data = load(out_file);

r2_gpfa_mean = mean(data.r2_gpfa, 2);

out_file = fullfile(results_dir, 'smoothed_neural_decode_r2only.mat');
data = load(out_file);

r2_smoothed_mean = mean(data.r2_smoothed, 2);

mat = cat(2, r2_smoothed_mean, r2_gpfa_mean, r2_single_mean, r2_stitched_mean);

%% Plot line plot with decode types on x-axis and R^2 on y axis, connected by dataset

clf;
figSize([4.5 5]);

h = plot(mat', '-', 'Color', grey(0.3));
TrialDataUtilities.Plotting.setLineOpacity(h, 0.7);
hold on;

x = repmat([1 2 3 4], size(mat, 1), 1);
h = scatter(x(:), mat(:), 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', grey(0.8), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.4);
% h = plot(mat', 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', grey(0.2), 'MarkerSize', 2);
% TrialDataUtilities.Plotting.setMarkerOpacity(h, 0.4);

h = plot(median(mat), 'r-', 'LineWidth', 2);
TrialDataUtilities.Plotting.setLineOpacity(h, 0.8); 


tickLabels = {sprintf('Smoothed\nNeural'), 'GPFA', ...
    sprintf('Single\nLFADS'), sprintf('Stitched\nLFADS')};

% single vs stitched bridge
e = 0.01;
yoff = 0.06;
yv = max(max(mat(:, [3 4]))) + yoff;
plot([3+e 4-e], [yv yv], 'k-', 'LineWidth', 1);
text(3.5, yv-0.01, '***', 'Background', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% single vs gpfa bridge
e = 0.01;
yv = max(max(mat(:, [2 3]))) + yoff + 0.05;
plot([2+e 3-e], [yv yv], 'k-', 'LineWidth', 1);
text(2.5, yv-0.01, '***', 'Background', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

k
e = 0.01;
yv = max(max(mat(:, [2 3]))) + 0.04;
plot([1+e 3-e], [yv yv], 'k-', 'LineWidth', 1);
text(2, yv-0.01, '***', 'Background', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% gpfa vs smoothed bridge
e = 0.01;
yv = max(max(mat(:, [1 2]))) + 0.04;
plot([1+e 2-e], [yv yv], 'k-', 'LineWidth', 1);
text(1.5, yv-0.01, '***', 'Background', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');


xlim([0.8 4.2]);
ylim([-0 1]);

hold off;
ylabel(sprintf('R^2 cross-validated\nkinematic predictions'));
aa = AutoAxis();
aa.addAutoAxisY();
aa.addTickBridge('x', 'tick', 1:4, 'tickLabel', tickLabels); 
aa.axisMarginTop = 0.2;
aa.axisMarginLeft = 1.3;
aa.axisMarginRight = 0.3;
aa.axisMarginBottom = 1;
aa.update();

%%

saveFigure(fullfile(fig_dir, 'kinematicCrossComparison'), 'painters', true, 'upsample', 4)

%% Summary stats

med = median(mat);
qtr = prctile(mat, [25 75], 1);

varNames = {'median', 'prc25', 'prc75'};
rowNames = {'smoothed', 'gpfa', 'single_lfads', 'multi_lfads'};

stats = array2table([med; qtr]', 'RowNames', rowNames', 'VariableNames', varNames);

diary(fullfile(results_dir, 'results_allDecoding.txt'));
stats

fid = fopen(fullfile(results_dir, 'results_kinematicSummary.txt'), 'w');

fprintf(fid,'Stitched vs single: p = %g sign rank, mean delta %g\n', ...
    signrank(r2_stitched_mean, r2_single_mean), mean(r2_stitched_mean - r2_single_mean));

fprintf(fid, 'Smoothed vs single: p = %g sign rank, mean delta %g\n', ...
    signrank(r2_smoothed_mean, r2_single_mean), mean(r2_single_mean - r2_smoothed_mean));
fprintf(fid, 'GPFA vs single: p = %g sign rank, mean delta %g\n', ...
    signrank(r2_gpfa_mean, r2_single_mean), mean(r2_single_mean - r2_gpfa_mean));


fprintf(fid, 'GPFA vs smoothed neural: p = %g sign rank, mean delta %g\n', ...
    signrank(r2_gpfa_mean, r2_smoothed_mean), mean(r2_gpfa_mean - r2_smoothed_mean));

fclose(fid);

diary off;

