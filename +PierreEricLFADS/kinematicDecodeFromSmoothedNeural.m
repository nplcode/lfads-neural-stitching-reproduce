function out = kinematicDecodeFromSmoothedNeural(rc, out_path, recompute)
    
    if nargin < 2
        out_path = fullfile(rc.path, 'results');
    end
    if nargin < 3
        recompute = false;
    end
    
     % check for existing
    r2_file = fullfile(results_dir, 'smoothed_neural_decode_r2only.mat');
    kinematics_file = fullfile(out_path, sprintf('kinematics_smoothed_neural.mat'));
    if ~recompute && exist(r2_file, 'file') && exist(kinematics_file, 'file')
        out = load(kinematics_file);
        d = load(r2_file);
        out.r2_smoothed = d.r2_smoothed;
        return;
    end
    
    ra = rc.findRuns('all');
    D = ra.nDatasets;
    
    %% Peri-move period
    
    opts = struct();
    opts.lag = 90;
    opts.align = 'GoCue';
    opts.tStart = 0; % was 200
    opts.binWidth = 20;
    
    
    opts_smooth = keepfields(opts, {'lag', 'tStart', 'align', 'binWidth'});
    opts_smooth.source = 'smoothed_neural';
    opts_smooth.neural_smooth = 40;
    
    Kfold = 5;
    
    %%
    % build the T struts using gpfa_xorth
    Tin_smoothed = ra.buildTStructs(opts_smooth);
    
    %% Do decoding one dataset at a time
    
    r2_smoothed = nan(ra.nDatasets, 2);
    
    prog = ProgressBar(ra.nDatasets, 'Decoding from smoothed neural');
    Tout_smoothed = cellvec(ra.nDatasets);
    
    for iD = 1:ra.nDatasets
        prog.update(iD);
        
        Tin_this = Tin_smoothed{iD};
        Tout_this = PierreEricLFADS.Decoding.kFoldXvalOLE(Tin_this, Kfold);
        r2_smoothed(iD, :) = Analysis.calcR2_TT(Tin_this, Tout_this)';
        Tout_smoothed{iD} = Tout_this;
    end
    prog.finish();
    
    
    %% Save final r2 values
    
    out_file = fullfile(out_path, 'smoothed_neural_decode_r2only.mat');
    save(out_file, 'r2_smoothed');
    
    %% Save decoded kinematics too
    
    if ~exist(out_path, 'dir')
        mkdirRecursive(out_path);
    end
    
    prog = ProgressBar(D, 'Saving decodes');
    for iDS = 1:D
        prog.update(iDS);
        
        out_file = fullfile(out_path, sprintf('kinematics_smoothed_neural_%s.mat', ra.datasets(iDS).name));
        T_smoothed_neural = Tout_smoothed{iDS}; %#ok<NASGU>
        save(out_file, 'T_smoothed_neural');
    end
    prog.finish();
    
    out_file = fullfile(out_path, 'kinematics_smoothed_neural.mat');
    save(out_file, 'Tin_smoothed', 'Tout_smoothed');
    
    %% assemble outputs
    
    out.Tin_smoothed = Tin_smoothed;
    out.Tout_smoothed = Tout_smoothed;
    out.r2_smoothed = r2_smoothed;

end
