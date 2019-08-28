function r2_gpfa = doDecodeFromGPFASpecificHPs(results_dir, rc, gpfa_binsize, gpfa_numlatents)

    if nargin < 2
        out_path = fullfile(rc.path, 'results');
    end
    if nargin < 3
        recompute = false;
    end
        
    % check for existing
    r2_file = fullfile(results_dir, 'gpfa_decode_r2only.mat');
    kinematics_file = fullfile(out_path, sprintf('kinematics_gpfa.mat'));
    if ~recompute && exist(r2_file, 'file') && exist(kinematics_file, 'file')
        out = load(kinematics_file);
        d = load(r2_file);
        out.r2_gpfa = d.r2_gpfa;
        return;
    end
    
    ra = rc.findRuns('all');
    D = ra.nDatasets;

%% run GPFA

ra.runGPFA(gpfa_binsize, gpfa_numlatents, false);

%% Peri-move period

opts = struct();
opts.lag = 90;
opts.align = 'GoCue';
opts.tStart = 0; % was 200
opts.binWidth = 20;

opts_gpfa = keepfields(opts, {'lag', 'tStart', 'align', 'binWidth'});
opts_gpfa.source = 'gpfa_xorth';
opts_gpfa.neuralBinSize = gpfa_binsize;

Kfold = 5;

% get sequence including gpfa results
ra.addGPFAResultsToTrialData();

% build the T struts using gpfa_xorth
Tin_gpfa = ra.buildTStructs(opts_gpfa);

% Tin_smoothed = ra.buildTStructs(opts_smooth);

%% Do decoding one dataset at a time

r2_gpfa = nan(ra.nDatasets, 2);

% prog = ProgressBar(ra.nDatasets, 'Decoding from smoothed neural and GPFA');
Tout_gpfa = cellvec(ra.nDatasets);

for iD = 1:ra.nDatasets
    % prog.update(iD);
    Tin_this = Tin_gpfa{iD};
    Tout_this = PierreEricLFADS.Decoding.kFoldXvalOLE(Tin_this, Kfold);
    r2_gpfa(iD, :) = Analysis.calcR2_TT(Tin_this, Tout_this)';
    Tout_gpfa{iD} = Tout_this;
end
% prog.finish();

desc = sprintf('bin%dms_latents%d', gpfa_binsize, gpfa_numlatents);

%% Save final r2 values

out_path = fullfile(results_dir, 'gpfa_sweep');
if ~exist(out_path, 'dir')
    mkdir(out_path);
end

out_file = fullfile(out_path, sprintf('gpfa_decode_r2only_%s.mat', desc));
save(out_file, 'r2_gpfa');

%% Save decoded kinematics too

% prog = ProgressBar(D, 'Saving decodes');
% for iDS = 1:D
%     % prog.update(iDS);
%     out_file = fullfile(out_path, sprintf('kinematics_gpfa_%s_%s.mat', desc, ra.datasets(iDS).name));
%     T_gpfa = Tout_gpfa{iDS};
%     save(out_file, 'T_gpfa');
% end
% prog.finish();

out_file = fullfile(out_path, sprintf('kinematics_gpfa_%s.mat', desc));
save(out_file, 'Tin_gpfa', 'Tout_gpfa');
