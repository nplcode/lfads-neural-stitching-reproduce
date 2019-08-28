function out = kinematicDecodeFromFactors(rc, results_dir, recompute, saveFigures)
    
    if nargin < 2
        results_dir = fullfile(rc.path, 'results');
    end
    if nargin < 3
        recompute = false;
    end
    if nargin < 4
        saveFigures = false;
    end
        
    % check for existing
    kinematics_stitched_file = fullfile(results_dir, 'kinematics_lfadsStitchedFromFactors.mat');
    kinematics_single_file = fullfile(results_dir, 'kinematics_lfadsSingleFromFactors.mat');
    r2_file = fullfile(results_dir, 'lfads_kinematics_fromFactors.mat');
    
    if ~recompute && exist(r2_file, 'file') && exist(kinematics_stitched_file, 'file') && exist(kinematics_single_file, 'file')
        out = load(kinematics_stitched_file);
        d = load(kinematics_single_file);
        out.Tin_single = d.Tin_single;
        out.Tout_single = d.Tout_single;
        d = load(r2_file);
        out.r2_single_factors = d.r2_single_factors;
        out.r2_stitched_factors = d.r2_stitched_factors;

    else
    
        ra = rc.findRuns('all');
        
        %% Peri-move period
        opts = struct();
        opts.lag = 90;
        opts.align = 'GoCue';
        opts.tStart = 0; % was 200
        opts.binWidth = 20;
        opts.source = 'factors';
        
        Kfold = 5;
        
        ra = rc.findRuns('all');
        
        optsT = keepfields(opts, {'lag', 'source', 'tStart', 'align', 'binWidth'});
        
        %% run cross-validated decoding on each single dataset rates individually
        
        %     lambdas = [0, 10.^(linspace(-1, 1, 11))];
        
        D = ra.nDatasets;
        prog = ProgressBar(D, 'Decoding datasets');
        
        r2_single_factors = deal(nan(D, 2));
        Tout_single = cellvec(D);
        %     lambda_best_single = nanvec(D);
        Tin_single = cellvec(D);
        for dsIdx = 1:D
            prog.update(dsIdx);
            r = rc.runs(dsIdx);
            Tsingle = r.buildTStructs(optsT);
            Tsingle = Tsingle{1};
            
            Tin_single{dsIdx} = Tsingle;
            
            %     [, Tout_single{dsIdx}] = ...1
            %         Analysis.Decoding.findBest_kFoldXvalOLE_l2reguralize(Tsingle, Kfold, 0);
            
            [Tout_single{dsIdx}] = ...
                PierreEricLFADS.Decoding.kFoldXvalOLE(Tsingle, Kfold);
            
            r2_single_factors(dsIdx, :) = Analysis.calcR2_TT(Tsingle, Tout_single{dsIdx})';
        end
        prog.finish();
        
        %% run single decoder on stitched datasets
        
        D = ra.nDatasets;
        Tin_stitched = ra.buildTStructs(optsT);
        
        %% Find cross validated decoder for all stitched datasets
        [Tstitched_cat, dsIdx] = TensorUtils.catWhich(1, Tin_stitched{:});
        Tstitched_cat = makecol(Tstitched_cat);
        
        % [r2_combined, Tstitched_catout_reg] = Analysis.Decoding.findBest_kFoldXvalOLE(Tstitched_cat, Kfold, 0);
        [Tstitched_catout, model] = PierreEricLFADS.Decoding.kFoldXvalOLE(Tstitched_cat, Kfold);
        
        % split back into individual days and compute per day R2
        
        Tout_stitched = TensorUtils.splitAlongDimensionByIndex(makecol(Tstitched_catout), 1, dsIdx);
        
        r2_stitched_factors = nan(D, 2);
        for iD = 1:D
            r2_stitched_factors(iD, :) = Analysis.calcR2_TT(Tin_stitched{iD}, Tout_stitched{iD})';
        end
        
        %% Save results
        
        save(r2_file, 'r2_single_factors', 'r2_stitched_factors');
        
        %% Save individual predictions
        
        out_path = fullfile(results_dir, 'kinematic_predictions');
        if ~exist(out_path, 'dir')
            mkdir(out_path);
        end
        
        prog = ProgressBar(D, 'Saving kinematics');
        for iDS = 1:D
            prog.update(iDS);
            out_file = fullfile(out_path, sprintf('kinematics_true_%s.mat', ra.datasets(iDS).name));
            T_true = Tin_single{iDS};
            save(out_file, 'T_true');
            
            out_file = fullfile(out_path, sprintf('kinematics_lfadsSingleFromFactors_%s.mat', ra.datasets(iDS).name));
            T_lfads_single = Tout_single{iDS};
            r2 = r2_single_factors(iDS, :);
            save(out_file, 'T_lfads_single', 'r2');
            
            out_file = fullfile(out_path, sprintf('kinematics_lfadsStitchedFromFactors_%s.mat', ra.datasets(iDS).name));
            T_lfads_stitched = Tout_stitched{iDS};
            r2 = r2_stitched_factors(iDS, :);
            save(out_file, 'T_lfads_stitched', 'r2');
        end
        prog.finish();
        
        debug('Saving predictions for all\n')
        
        save(kinematics_stitched_file, 'Tin_stitched', 'Tout_stitched');
        
        save(kinematics_single_file, 'Tin_single', 'Tout_single');
        
        out.Tin_stitched = Tin_stitched;
        out.Tout_stitched = Tout_stitched;
        out.Tin_single = Tin_single;
        out.Tout_single = Tout_single;
        out.r2_single_factors = d.r2_single_factors;
        out.r2_stitched_factors = d.r2_stitched_factors;
    end

    %% Plot R2 comparison
    
    figure(1); clf
    hx = scatter(r2_single_factors(:, 1), r2_stitched_factors(:, 1), 12, 'filled', 'MarkerEdgeColor', 'w', 'MarkerEdgeAlpha', 0.4, 'LineWidth', 0.1);
    hold on
    hy = scatter(r2_single_factors(:, 2), r2_stitched_factors(:, 2), 12, 'filled', 'MarkerEdgeColor', 'w', 'MarkerEdgeAlpha', 0.4, 'LineWidth', 0.1);
    TrialDataUtilities.Plotting.showInLegend(hx, 'X');
    TrialDataUtilities.Plotting.showInLegend(hy, 'Y');
    % h = legend(gca, 'show'); h.Location = 'SouthEast';
    axis equal;
    set(gca, 'TickDir', 'out');
    title('Kinematic Predictions');
    % title(sprintf('Comparison lag %d bin %d source %s', opts.lag, opts.binWidth, opts.source), 'Background', 'none');
    set(xlabel('Single-session R^2'), 'BackgroundColor', 'none');
    set(ylabel('Multi-session R^2'), 'BackgroundColor', 'none');
    xlim([-0.2 1]);
    ylim([-0.2 1]);
    hold on;
    hl = TrialDataUtilities.Plotting.identityLine;
    uistack(hl, 'bottom');
    % niceGrid;
    
    figSize([5 5]);
    ax = AutoAxis.replace();
    ax.addColoredLabels({'x', 'y'}, [hx.CData; hy.CData], 'posY', AutoAxis.PositionType.Bottom, 'posX', AutoAxis.PositionType.Right);
    ax.update();
    hold off;

    if saveFigures
        saveFigure(fullfile(results_dir, 'kinematics_r2_compare_single_vs_stitching_fromFactors'), 'painters', true, 'upsample', 4);
    end
    
    %% Print sign rank p values
    
    fname = fullfile(results_dir, 'results_kinematicsFromFactors.txt');
    fid = fopen(fname, 'w');
    
    tee(fid, 'p value signrank x: %.8g\n', signrank(r2_stitched_factors(:, 1), r2_single_factors(:, 1)));
    tee(fid, 'p value signrank y: %.8g\n', signrank(r2_stitched_factors(:, 2), r2_single_factors(:, 2)));
    
    tee(fid, 'median delta %.8g\n', median(r2_stitched_factors(:) - r2_single_factors(:)));
    
    tee(fid, 'delta of medians %.8g\n', median(r2_stitched_factors(:)) - median(r2_single_factors(:)));
    
    fclose(fid);
    
    dbtype(fname);
end

function tee(fid, varargin)
    fprintf(varargin{:});
    fprintf(fid, varargin{:});
end
    
