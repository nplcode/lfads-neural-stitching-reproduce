function out = kinematicDecodeFromGPFA(rc, out_path, recompute)

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
    
    %% Load best GPFA params from sweep
    
    gpfaSweep_file = fullfile(results_dir, 'gpfa_sweep/gpfa_sweep_r2only.mat');
    d = load(gpfaSweep_file);
    gpfa_binSize = d.best_binSize;
    gpfa_numlatents = d.best_numLatents;
    
    %% run GPFA
    
    % gpfa_binsize = 5;
    % gpfa_numlatents = 20;
    ra.runGPFA(gpfa_binsize, gpfa_numlatents, false);
    
    %% Peri-move period
    
    opts = struct();
    opts.lag = 90;
    opts.align = 'GoCue';
    opts.tStart = 0; % was 200
    opts.binWidth = 20; % this is what we will sample the kinematics and the gpfa signal at
    
    opts_gpfa = keepfields(opts, {'lag', 'tStart', 'align', 'binWidth'});
    opts_gpfa.source = 'gpfa_xorth';
    opts_gpfa.neuralBinSize = gpfa_binsize;
    
    Kfold = 5;
    
    %% Build GPFA T structs
    
    % get sequence including gpfa results
    ra.addGPFAResultsToTrialData();
    
    % build the T struts using gpfa_xorth
    Tin_gpfa = ra.buildTStructs(opts_gpfa);
    
    %% Do decoding one dataset at a time
    
    r2_gpfa = nan(ra.nDatasets, 2);
    
    prog = ProgressBar(ra.nDatasets, 'Decoding from smoothed neural and GPFA');
    Tout_gpfa = cellvec(ra.nDatasets);
    
    for iD = 1:ra.nDatasets
        prog.update(iD);
        Tin_this = Tin_gpfa{iD};
        Tout_this = PierreEricLFADS.Decoding.kFoldXvalOLE(Tin_this, Kfold);
        r2_gpfa(iD, :) = Analysis.calcR2_TT(Tin_this, Tout_this)';
        Tout_gpfa{iD} = Tout_this;
    end
    prog.finish();
    
    %% Save final r2 values
    save(r2_file, 'r2_gpfa');
    
    %% Save decoded kinematics too
    
    out_path = fullfile(results_dir, 'kinematic_predictions');
    if ~exist(out_path, 'dir')
        mkdir(out_path);
    end
    
    prog = ProgressBar(D, 'Saving decodes');
    for iDS = 1:D
        prog.update(iDS);
        out_file = fullfile(results_dir, sprintf('kinematics_gpfa_%s.mat', ra.datasets(iDS).name));
        T_gpfa = Tout_gpfa{iDS};
        save(out_file, 'T_gpfa');
    end
    prog.finish();
    
    save(kinematics_file, 'Tin_gpfa', 'Tout_gpfa');

    out.Tin_gpfa = Tin_gpfa;
    out.Tout_gpfa = Tout_gpfa;
    out. 