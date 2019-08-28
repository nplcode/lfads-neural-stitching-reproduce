classdef Run < LFADS.Run
    properties
        trialDataSet % loaded from database
        gpfaSequenceData % loaded from GPFA run
    end
    
    properties(Dependent)
        pathGPFAOutput
    end
    
    methods % needed for LFADS prep
        function r = Run(varargin) 
           r@LFADS.Run(varargin{:});
        end
        
        function seq = convertDatasetToSequenceStruct(r, dataset, mode)
            data = dataset.loadData();
            data = data.data;
            
            preKeep = double(r.params.preKeep);
            postKeep = double(r.params.postKeep);
            
            runID = [r.name '__' r.paramsString '__' dataset.name];
            
            if strncmpi(r.datasetCollection.context, 'Eric', 4)
                condField = 'targetDirectionName';
            else
                condField = 'conditionDesc';
            end
            
            seq = [];
            for it = 1:data.nTrials
                % don't keep the ultra-short trials that reach target early
                %                 if postGoTimeToKeep > data.TargetAcquired(it)
                %                     continue;
                %                 end
                % some GoCue's are nan...
                if isnan(data.(r.params.align)(it))
                    error('GoCue is nan');
                end
                
                seq(it).runID = runID; %#ok<*AGROW>
                seq(it).trialID = data.trialId(it);
                seq(it).saveTag = dataset.saveTags;
                seq(it).conditionId = data.(condField){it};
                
                if strcmpi(r.datasetCollection.context, 'Eric')
                    seq(it).peakSpeed2 = data.PeakSpeed2(it);
                    seq(it).peakSpeed3 = data.PeakSpeed3(it);
                    seq(it).rt = data.RT(it);
                    seq(it).delay = data.GoCue(it) - data.TargetOnset(it);
                end
                
                % data is aligned to target onset
                % we'll just take time starting at target onset
                timeIndsToKeep = data.(r.params.align)(it) + (-preKeep:postKeep-1);
                spStart = argmin(abs(data.spikeRasters_time - timeIndsToKeep(1)));
                spikeRasterIndsToKeep = spStart + (1:numel(timeIndsToKeep));
%                 [y_time, spikeRasterIndsToKeep] = intersect(data.spikeRasters_time, ...
%                     timeIndsToKeep);
                seq(it).y = squeeze(data.spikeRasters(it, spikeRasterIndsToKeep, ...
                    :))';
                if any(isnan(seq(it).y(:)))
%                     warning('nans in the y matrix');
                    seq(it).y(isnan(seq(it).y)) = 0;
                end
                seq(it).y_time = -preKeep:postKeep-1;
                seq(it).T = size(seq(it).y,2);
                if seq(it).T == 0
                    error('Issue with data in spike raster')
                end
                
                % store down info about target name
                seq(it).targetDirectionName = data.targetDirectionName{it};
                
                if strcmpi(r.datasetCollection.context, 'Dan')
                    seq(it).conditionDesc = data.conditionDesc{it};
                end
                
                % store down some hand kinematics
                x=horzcat(data.handKinematics{it,:});
                [hand_time, handKinematicsIndsToKeep] = intersect(data.handKinematics_time{it}, ...
                    timeIndsToKeep);
                seq(it).hand_time = double(hand_time - data.(r.params.align)(it));
                seq(it).handKinematics = x(handKinematicsIndsToKeep,:)';
                
                seq(it).binWidthMs = 1;
                seq(it).params.dtMS = 1;
                seq(it).params.runID = runID;
                seq(it).subject = dataset.subject;
                seq(it).date = dataset.datenum;
            end
            
            nTrialsKeep = r.params.nTrialsKeep;
            if numel(seq) > nTrialsKeep
                seq = seq(1:nTrialsKeep);
            end
            
            if r.params.minDelay > 0
                delay = cat(1, seq.delay);
                seq = seq(delay >= r.params.minDelay);
            end
            
            if r.params.maxDelay < Inf
                delay = cat(1, seq.delay);
                seq = seq(delay <= r.params.maxDelay);
            end
            
            % check for no spike units
            ycat = cat(3, seq.y); % nCh x T x nTrials
            nSpikesTotal = sum(sum(ycat, 2), 3);
            if any(nSpikesTotal == 0)
                debug('Removing %d units with no spikes\n', nnz(nSpikesTotal == 0));
                mask = nSpikesTotal > 0;
                for it = 1:numel(seq)
                    seq(it).y = seq(it).y(mask, :);
                    seq(it).unitMask = mask;
                end
            end  
            
            seq = LFADS.Utils.makecol(seq);
        end
    end
    
    methods % Working with LFADS generated data with TrialData
        
        function tdSet = loadTrialDataFromDatasetCollection(r)
            % pierreSettingsEric
            % db = LFADS_PierreExport.Database.loadDb();
            
            if isempty(r.datasetCollection.trialDataSet)
                error('Call .loadTrialDataFromDatabase on datasetCollection');
            end
            tdSetFull = r.datasetCollection.trialDataSet;
            tdSet = cell(r.nDatasets, 1);
            dsidx = r.datasetIndsInCollection;
            
            prog = ProgressBar(r.nDatasets, 'Preparing trialData from datasetCollection');
            for i = 1:r.nDatasets
                prog.update(i);
               
                td = tdSetFull{dsidx(i)};
                if ~isempty(td)
                    td = TrialDataConditionAlign(td);
                    
                    delay = td.getParam('delay');
                    mask = delay >= r.params.minDelay & delay <= r.params.maxDelay;
                    td = td.selectTrials(mask);
                     
                    % start with same alignment of neural data
                    td = td.start(r.params.align, -r.params.preKeep).stop(r.params.align, r.params.postKeep);
                    td = td.groupBy('targetDirectionName').setAttributeDisplayAs('targetDirectionName', 'Target');
                    tdSet{i} = td;
                else
                    error('Could not locate run %d', i);
                end
            end
            prog.finish();
            
            r.trialDataSet = tdSet;
        end
        
        function tdSet = addPosteriorMeansToTrialData(r)
            if isempty(r.trialDataSet)
                r.loadTrialDataFromDatasetCollection();
            end
            r.loadPosteriorMeans();
            
            timeField = 'posteriorMeans_time';
            
            prog = ProgressBar(r.nDatasets, 'Merging posterior mean data into trialData for each dataset');
            tdSet = cellvec(r.nDatasets);
            for i = 1:r.nDatasets
                prog.update(i);
                td =  r.trialDataSet{i};
                pm = r.posteriorMeans(i);
                
                % data must be nTrials x nTime x nChannels tensor
                td = td.reset.zero(r.params.align);
                td = td.dropAnalogChannelGroup({'controllerOutputs', 'factors', 'generatorStates', 'rates'});
                td = td.dropChannel('generatorIC');
                
                td = td.addAnalogChannelGroup('factors', genNames('f', pm.nFactors), ...
                    permute(pm.factors, [3 2 1]), pm.time, 'timeField', timeField, 'isAligned', true);
                td = td.addAnalogChannelGroup('generatorStates', genNames('g', pm.nGeneratorUnits), ...
                    permute(pm.generator_states, [3 2 1]), pm.time, 'timeField', timeField, 'isAligned', true);
                td = td.addAnalogChannelGroup('rates', genNames('r', pm.nNeurons), ...
                    permute(pm.rates, [3 2 1]), pm.time, 'timeField', timeField, 'isAligned', true);
                
                if pm.nControllerOutputs > 0
                    td = td.addAnalogChannelGroup('controllerOutputs', genNames('co', pm.nControllerOutputs), ...
                        permute(pm.controller_outputs, [3 2 1]), pm.time, 'timeField', timeField, 'isAligned', true);
                end
                
                td = td.addVectorParamAccessAsMatrix('generatorIC', TensorUtils.splitAlongDimension(pm.generator_ics, 2)');
                
                tdSet{i} = td;
            end
            
            % automatically align trial data the same way as the LFADS data was
            % prepared
            tdSet = r.setupTrialDataAsLFADS(tdSet);
            r.trialDataSet = tdSet;
            
            function names = genNames(pre, n)
                names = arrayfun(@(i) sprintf('%s%i', pre, i), (1:n)', 'UniformOutput', false);
            end
        end
        
        function tdSet = getTrialDataWithVirtualPopulationRate(r, dsIdx, W, b)
            timeField = 'posteriorMeans_time';
            
            % concatenated readout matrices for each dataset
            ro = r.loadReadoutMatricesByDataset();
            
            if ~exist('W', 'var')
                W = cat(1, ro.rates_W); % N_all x Far')
                b = cat(1, ro.rates_b); % N_all x 1
            end
            
            % row normalize as is done in lfads code
            Wnorm = W ./ sqrt(sum(W.^2, 2));
            N_all = size(Wnorm, 1);
            
            vrNames = genNames('vr', N_all);
            dsIdx = TensorUtils.vectorMaskToIndices(dsIdx);
            
            prog = ProgressBar(numel(dsIdx), 'Adding virtual population rates to trialData for each dataset');
            tdSet = cellvec(numel(dsIdx));
            for i = 1:numel(tdSet)
                idx = dsIdx(i);
                prog.update(idx);
                td =  r.trialDataSet{idx};
                pm = r.posteriorMeans(idx);
                
                % data must be nTrials x nTime x nChannels tensor
                td = td.reset.zero(r.params.align);
                td = td.dropAnalogChannelGroup({'virtualRates'});
                
                virtualRates = exp(TensorUtils.linearCombinationAlongDimension(pm.factors, 1, Wnorm));
                
                td = td.addAnalogChannelGroup('virtualRates', vrNames, ...
                    permute(virtualRates, [3 2 1]), pm.time, 'timeField', timeField, 'isAligned', true);
                tdSet{i} = td;
            end
            
            % automatically align trial data the same way as the LFADS data was
            % prepared
            tdSet = r.setupTrialDataAsLFADS(tdSet);
            
            function names = genNames(pre, n)
                names = arrayfun(@(i) sprintf('%s%i', pre, i), (1:n)', 'UniformOutput', false);
            end
        end
        
        function tdSet = setupTrialDataAsLFADS(r, tdSet)
            for i = 1:numel(tdSet)
                td =  tdSet{i};
                td = td.start(r.params.align, -r.params.preKeep).stop(r.params.align, r.params.postKeep);
                td = LFADS_PierreExport.Condition.groupCenterOut(td);
                td = td.setAttributeDisplayAs('targetDirectionName', 'Target');
                tdSet{i} = td;
            end
        end
        
        function Tset = buildTStructs(r, varargin)
            % opts can specify :
            %   opts.lag - how much to lag the neural and kinematic data
            %     (each trial will be trimmed by lag milliseconds)
            % tStart (scalar) : where to start relative to r.params.align,
            % if NaN, use the first output from LFADS
            %   opts.neuralBinSizeMS - what is the neural data currently binned at?
            %   opts.neuralFieldName - specify the field name of the neural
            %     data. Default is 'y'

            p = inputParser();
            p.addParameter('datasetMask',truevec(r.nDatasets), @(x) true);
            p.addParameter('source', 'neural', @(x) ismember(x, {'neural', 'smoothed_neural', 'rates', 'virtualRates', ...
                'factors', 'generatorStates', 'gpfa_xorth'}));
            p.addParameter('neural_smooth', 40, @isscalar); % for smoothed_neural, SD of the Gaussian
            p.addParameter('rates_W', [], @ismatrix); % for virtualRates, defaults to full set of virtual neurons
            p.addParameter('rates_b', [], @isvector); % for virtualRates, defaults to full set of virtual neurons
            p.addParameter('align', 'GoCue', @ischar);
            p.addParameter('binWidth', 20, @isscalar);
            p.addParameter('tStart', NaN, @isscalar);
            p.addParameter('tStop', NaN, @isscalar);
            p.addParameter('lag', 0, @isscalar);
            p.addParameter('neuralFieldName', 'y', @ischar);
            p.addParameter('neuralBinSize', 1, @isscalar);
            p.parse(varargin{:});
            
            tLag = p.Results.lag;
            binWidth = p.Results.binWidth;
            source = p.Results.source;
            align = p.Results.align;
            
            dsIdx = TensorUtils.vectorMaskToIndices(p.Results.datasetMask);
            
            % take the limits of the posterior means in the appropriate
            % alignment
            td = r.trialDataSet{dsIdx(1)};
   
            prog = ProgressBar(r.nDatasets, 'Building T structs');
            for iDS = 1:numel(dsIdx)
                prog.update(iDS);
                td = r.trialDataSet{dsIdx(iDS)};
                
                td = td.unalign.zero(align);
            
                [data, tvec] = td.getAnalogChannelGroupAsTensor('rates', 'minTrialFraction', 1);
                tStart = nanmax(p.Results.tStart, tvec(1));
                tStop = nanmin(p.Results.tStop, tvec(end));
                
                % round to nearest binWidth multiple
                tStart = ceil(tStart / binWidth) * binWidth;
                tStop = floor(tStop / binWidth) * binWidth;
                
                % align trial data
                td = td.unalign.start(align, tStart).stop(align, tStop);
                
                % lag kinematics relative to neural data 
                % Trials x Time x 4 channels --> 4 x Time x Trials
                [kinData, time] = td.lag(tLag).getAnalogMultiAsTensor(...
                    {'handX', 'handY', 'handVelocityX', 'handVelocityY'}, 'timeDelta', binWidth);
                if any(isnan(kinData(:)))
                    warning('%d nans in kinematics data', nnz(isnan(kinData(:))));
                end
                kinematics = permute(kinData, [3 2 1]);
                
                % append ones as lasst channel
                kinematics = TensorUtils.expandAlongDims(kinematics, 1, 1, 1);
                
                % split by trials
                kinematics = squeeze(TensorUtils.splitAlongDimension(kinematics, 3));
                
                switch source
                    case 'neural'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        decode = permute(td.getSpikeBinnedCounts(td.listSpikeChannels, 'binWidthMs', binWidth), [3 2 1]);
                    
                    case 'smoothed_neural'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        sf = GaussianSpikeFilter('sigma', p.Results.neural_smooth);
                        sf.timeDelta = binWidth;
                        decode = permute(td.getSpikeRateFilteredAsMatrix(td.listSpikeChannels, 'spikeFilter', sf), [3 2 1]);
                        mask = isnan(decode);
                        if any(mask(:))
                            warning('Setting %d values from NaN to 0', nnz(mask));
                            decode(mask) = 0;
                        end
                            
                    case 'rates'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        decode = permute(td.getAnalogChannelGroupAsTensor('rates', 'timeDelta', binWidth), [3 2 1]);
                    
                    case 'virtualRates'
                        % Trials x Time x NeuronsAll --> Neurons x Time x Trials
                        factorTensor = permute(td.getAnalogChannelGroupAsTensor('factors', 'timeDelta', binWidth), [3 2 1]);
                        
                        % concatenated readout matrices for each dataset
                        ro = r.loadReadoutMatricesByDataset();
                        if isempty(p.Results.rates_W)
                            W = cat(1, ro.rates_W); % N_all x Factors
                            b = cat(1, ro.rates_b); % N_all x 1
                        else
                            W = p.Results.rates_W;
                            b = p.Results.rates_b;
                        end
                        Wnorm = W ./ sqrt(sum(W.^2, 2)); % row normalize as is done in lfads code  
            
                        % rates = exp(W*factors + b)
                        decode = exp(TensorUtils.linearCombinationAlongDimension(factorTensor, 1, Wnorm) + b);
                    
                    case 'factors'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        decode = permute(td.getAnalogChannelGroupAsTensor('factors', 'timeDelta', binWidth), [3 2 1]);
                    
                    case 'generatorStates'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        decode = permute(td.getAnalogChannelGroupAsTensor('generatorStates', 'timeDelta', binWidth), [3 2 1]);
                    
                    case 'gpfa_xorth'
                        % Trials x Time x nGPFA --> nGPFA x Time x Trials
                        decode = permute(td.getAnalogChannelGroupAsTensor('gpfa_xorth', 'timeDelta', binWidth), [3 2 1]);
                end
                
                % split by trials
                decode = squeeze(TensorUtils.splitAlongDimension(decode, 3));
                
                % number of timepoints
                nTimepoints = cellfun(@(x) size(x, 2), decode, 'UniformOutput', false);
                
                valid = num2cell(td.valid);
                
                conditions = num2cell(td.conditionIdx);
                
                % split by trials and assign into struct
                Tset{iDS} = struct('valid', valid, 'X', kinematics, 'Z', decode, 'time', time, 'T', nTimepoints, 'dt', p.Results.binWidth, 'condition', conditions);
            end
            
            prog.finish();
        end
        
        function tdSet = addDecodeTStructToTrialData(r, Tset, channelPrefix, varargin)
            % uses field xk to produce x and y velocity channels
            if isempty(r.trialDataSet)
                r.loadTrialDataFromDatasetCollection();
            end
            
            if ~iscell(Tset)
                assert(r.nDatasets == 1);
                Tset = {Tset};
            end
            
            timeField = sprintf('%s_time', channelPrefix);
            pxField = sprintf('%s_posx', channelPrefix);
            pyField = sprintf('%s_posy', channelPrefix);
            vxField = sprintf('%s_velx', channelPrefix);
            vyField = sprintf('%s_vely', channelPrefix);
            
            prog = ProgressBar(r.nDatasets, 'Merging Tstructs into trialData for each dataset');
            tdSet = cellvec(r.nDatasets);
            for i = 1:r.nDatasets
                prog.update(i);
                td =  r.trialDataSet{i};
                T = Tset{i};
                
                % data must be nTrials x 1 cell of nTime x nChannels
                % matrices
                dataThis = arrayfun(@(s) s.xk(1:4, :)', T, 'UniformOutput', false);
                timeThis = arrayfun(@(s) s.time, T, 'UniformOutput', false);
                td = td.reset.zero(r.params.align);
                td = td.dropAnalogChannelGroup(channelPrefix);
                
                td = td.addAnalogChannelGroup(channelPrefix, {pxField, pyField, vxField, vyField}, ...
                    dataThis, timeThis, 'timeField', timeField, 'isAligned', true);

                tdSet{i} = td;
            end
            
            % automatically align trial data the same way as the LFADS data was
            % prepared
            tdSet = r.setupTrialDataAsLFADS(tdSet);
            r.trialDataSet = tdSet;
        end
    end
    
    methods % GPFA smoothing of neural data
        function p = get.pathGPFAOutput(r)
            if isempty(r.runCollection)
                p = '';
            else
                p = fullfile(r.path, 'gpfaOutput');
            end
        end

        function p = getGpfaResultsFile(r, dtMS, nLatents)
            gpfaResultsDir = r.getGpfaResultsDir(dtMS, nLatents);
            p = fullfile(gpfaResultsDir, 'all_results.mat');
        end

        function gpfaResultsDir = getGpfaResultsDir(r, dtMS, nLatents)
            gpfaResultsDir = fullfile(r.pathGPFAOutput, ...
                                      sprintf('nLat_%03i_binSizeMS_%03i', nLatents, dtMS));
        end
        
        function resultsOut = loadGPFA(r, dtMS, nLatents)
            gpfaOutputFile = r.getGpfaResultsFile(dtMS, nLatents);
            tmp = load(gpfaOutputFile);
            resultsOut = tmp.resultsOut;
            r.gpfaSequenceData = resultsOut;
        end
        
        function resultsOut = runGPFA(r, dtMS, nLatents, deleteExistingResults)
        % function seq = doGPFA(r, dtMS, nLatents, deleteExistingResults)

            if ~exist('deleteExistingResults', 'var')
                deleteExistingResults = false;
            end
            
            % if sequence data is not yet loaded, do so now
            if isempty(r.sequenceData)
                r.loadSequenceData();
            end
            seqs = r.sequenceData;
            
            % need the training and validation inds
            if isempty(r.inputInfo)
                r.loadInputInfo();
            end
            
            trainInds = {r.inputInfo.trainInds}';
            validInds = {r.inputInfo.validInds}';

            xDim = nLatents;

            gpfaOutputFile = r.getGpfaResultsFile(dtMS, nLatents);

            % delete existing file if requested
            if deleteExistingResults
                if exist(gpfaOutputFile, 'file')
                    warning('Deleting existing GPFA results file %s', gpfaOutputFile);
                    delete(gpfaOutputFile)
                end
            end

            % load existing results if they exist
            if exist(gpfaOutputFile, 'file')
                fprintf('Loading previously-run GPFA results (this may take some time)\n');
                resultsOut = r.loadGPFA(dtMS, nLatents);
            else
                resultsOut = cell(numel(seqs), 1);
                binWarned = false;
                prog = ProgressBar(numel(seqs), 'Running GPFA on each dataset');
                for nseq = 1:numel(seqs)
                    prog.update(nseq);
                    seq = seqs{nseq};

                    rebin = dtMS / seq(1).params.dtMS;
                    if rebin ~= dtMS && ~binWarned
                        fprintf(['doGPFA: FYI, sequence data is binned ' ...
                                      'at %i ms, adjusting GPFA binning ' ...
                                      'appropriately, but this is ' ...
                                      'probably dumb\n'], seq(1).dtMS);
                        binWarned = true;
                    end


                    % GPFA is expecting data to be in a field
                    % called 'spikes'
                    [seq.spikes] = seq.y;
                    seq = rmfield(seq, 'y');

                    % GPFA is expecting each trial to be numbered
                    % with a "trialId"
                    tid = num2cell((1:numel(seq))');
                    [seq.trialId] = tid{:};

                    tic;
                    % run GPFA with a modified version that allows you to
                    % specify an output directory

                    thisSetGpfaResultsDir = fullfile(r.getGpfaResultsDir(dtMS, nLatents), ...
                                                     sprintf('set%03i', ...
                                                             nseq));
                    results = GPFA.neuralTraj_mod(1, seq, 'trainTrialIdx', trainInds{nseq}, ...
                                             'testTrialIdx', validInds{nseq}, 'saveDir', thisSetGpfaResultsDir,...
                                             'xDim', xDim, 'binWidth', rebin);
                    % postprocess - orthonormalization and "cleanup" (?)
                    [estParams2, seqTrain2, seqTest2] = GPFA.Util.postprocess(results);

                    % trim some of the unnecessary fields of the seqTrains
                    trains = {'seqTrain2','seqTest2'};
                    f2keep = {'trialId','T','y','xsm','xorth'};
                    for nn = 1:numel(trains)
                        tr = eval(trains{nn});
                        fs = fields(tr);
                        f2remove = setdiff(fs, f2keep);
                        for nf = 1:numel(f2remove)
                            tr = rmfield(tr, f2remove{nf});
                        end
                        eval(sprintf('%s = tr;', trains{nn}));
                    end

                    results.seqTrain = seqTrain2;
                    results.seqTest = seqTest2;
                    results.train_inds = trainInds{nseq};
                    results.valid_inds = validInds{nseq};
                    results.estParams = estParams2;
                    toc;

                    resultsOut{nseq} = results;
                end

                % check for necessity of v7.3
                varinfo=whos('resultsOut');
                saveopt='';
                if varinfo.bytes >= 2^31
                    saveopt='-v7.3';
                end

                save(gpfaOutputFile, 'resultsOut', ...
                     saveopt);
                prog.finish();
            end
           

            resultsOut = makecol(resultsOut);
            r.gpfaSequenceData = resultsOut;
        end
        
        function seqs = addGPFAResultsToSeq(r)
            % function seqs = addGPFAResultsToSeq(r)
            % returns a sequence that has posterior mean
            % values integrated

            if isempty(r.gpfaSequenceData)
                error(['Run.addGPFAResultsToSeq: first load results ' ...
                       'using r.doGPFA( ... )']);
            end

            if isempty(r.sequenceData) || numel(r.sequenceData) == 0
                r.sequenceData = r.loadSequenceData();
            end
            seqs = r.sequenceData;

            % iterate over datasets
            for iDS = 1:numel(seqs)
                gs = r.gpfaSequenceData{iDS};

                % training data
                for itr = 1:numel(gs.seqTrain)
                    ntr = gs.seqTrain(itr).trialId;
                    seqs{iDS}(ntr).gpfa_xorth = gs.seqTrain(itr).xorth;
                end

                % test data
                for itr = 1:numel(gs.seqTest)
                    ntr = gs.seqTest(itr).trialId;
                    seqs{iDS}(ntr).gpfa_xorth = gs.seqTest(itr).xorth;
                end
            end
        end
        
        function tdSet = addGPFAResultsToTrialData(r)
            % function seqs = addGPFAResultsToSeq(r)
            % returns a sequence that has posterior mean
            % values integrated

            if isempty(r.gpfaSequenceData)
                error(['first load results ' ...
                       'using r.doGPFA( ... )']);
            end

            if isempty(r.trialDataSet)
                r.loadTrialDataFromDatasetCollection();
            end
            tdSet = r.trialDataSet; 
            
            r.loadPosteriorMeans();

            % iterate over datasets
            prog = ProgressBar(r.nDatasets, 'Adding GPFA results to TrialData');
            for iDS = 1:r.nDatasets
                prog.update(iDS);
                gs = r.gpfaSequenceData{iDS};
                
                td = tdSet{iDS};
                pm = r.posteriorMeans(iDS);
                
                % data must be nTrials x nTime x nChannels tensor
                nGP = size(gs.seqTrain(1).xorth, 1);
                nTime = size(gs.seqTrain(1).xorth, 2);
                gpfa_xorth = nan(td.nTrials, nTime, nGP);
                
                binMs = gs.binWidth;
                preKeep = double(r.params.preKeep);
                postKeep = double(r.params.postKeep);
                gpfa_time = (-preKeep:binMs:postKeep-binMs)';
                assert(numel(gpfa_time) == nTime);

                % training data
                for itr = 1:numel(gs.seqTrain)
                    ntr = gs.seqTrain(itr).trialId;
                    gpfa_xorth(ntr, :, :) = gs.seqTrain(itr).xorth';
                end

                % test data
                for itr = 1:numel(gs.seqTest)
                    ntr = gs.seqTest(itr).trialId;
                    gpfa_xorth(ntr, :, :) = gs.seqTest(itr).xorth';
                end

                td = td.reset.zero(r.params.align);
                td = td.dropAnalogChannelGroup({'gpfa_xorth'}); 
                td = td.addAnalogChannelGroup('gpfa_xorth', genNames('gpfa', nGP), ...
                    gpfa_xorth, gpfa_time, 'timeField', 'gpfa_time', 'isAligned', true);
                
                tdSet{iDS} = td;
            end
            prog.finish();
            
            tdSet = r.setupTrialDataAsLFADS(tdSet);
            r.trialDataSet = tdSet;
            
            function names = genNames(pre, n)
                names = arrayfun(@(i) sprintf('%s%i', pre, i), (1:n)', 'UniformOutput', false);
            end
        end
    end
end
