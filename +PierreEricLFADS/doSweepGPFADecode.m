function out = doSweepGPFADecode(rc, out_path)

    ra = rc.findRuns('all');
    D = ra.nDatasets;
    

binSizes = [5 10 20]';
latentValues = (5:5:20)';
[B, L] = ndgrid(binSizes, latentValues);

N = numel(B(:));

prog = ProgressBar(N, 'Sweeping GPFA decodes');
% prog.enableParallel();

r2_cell = cell(size(B));

for i = 1:N
   gpfa_binsize = B(i);
   gpfa_numlatents = L(i);
   fname = ra.getGpfaResultsFile(gpfa_binsize, gpfa_numlatents);
   if ~exist(fname, 'file')
       ra.runGPFA(gpfa_binsize, gpfa_numlatents, false);
   end
    r2_cell{i} = PierreEricLFADS.Paper.doDecodeFromGPFASpecificHPs(results_dir, ra, B(i), L(i));
    prog.update(i);
end

prog.finish();

%% Pick the best GPFA params
r2_means = cellfun(@(x) mean(x(:)), r2_cell);

[~, idx] = max(r2_means(:));
best_binSize = B(idx);
best_numLatents = L(idx);


%% Save swept r2_values

out_path = fullfile(results_dir, 'gpfa_sweep');
if ~exist(out_path, 'dir')
    mkdir(out_path);
end

out_file = fullfile(out_path, 'gpfa_sweep_r2only.mat');
save(out_file, 'r2_cell', 'r2_means', 'binSizes', 'latentValues', 'B', 'L', 'best_binSize', 'best_numLatents');


%% Load the out_file

ld = load(out_file);
r2_cell = ld.r2_cell;
binSizes = ld.binSizes;
latentValues = ld.latentValues;

%% Plot performance


figure(1); clf;
for iB = 1:numel(binSizes)
    h = plot(latentValues, r2_means(iB, :), 'o-');
    hold on;
    TrialDataUtilities.Plotting.showInLegend(h, sprintf('%d ms bins', binSizes(iB)));
end

xlabel('Latent dimensionality');
ylabel('R^2');
ylim([0 1]);
xlim([min(binSizes) max(binSizes)]);

h = legend(gca, 'show');
legend boxoff;

AutoAxis.replace();
