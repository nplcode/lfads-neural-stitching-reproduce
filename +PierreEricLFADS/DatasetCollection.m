classdef DatasetCollection < LFADS.DatasetCollection
    properties
        context = 'Eric';
        trialDataSet % loaded from database
    end
    
    methods
        function ds = DatasetCollection(path)
            ds = ds@LFADS.DatasetCollection(path);
        end

        function autoDetectDatasets(dc)
            dc.clearDatasets;

            % automatically find all files within and build datasets
            files = dir(dc.path);
            for iF = 1:numel(files)
                if strncmp(files(iF).name, '.', 1), continue, end
                info = files(iF);
                [~, ~, ext] = fileparts(info.name);
                if ~strcmp(ext, '.mat'), continue; end
                ds = PierreEricLFADS.Dataset(dc, info.name);
            end
        end

        function filterHasHighSNRChannels(dc)
            dc.loadInfo();
            mask = arrayfun(@(ds) ds.nChannelsHighSNR > 0, dc.datasets);
            dc.filterDatasets(mask);
        end

        function filterBestSaveTagEachDate(dc)
            % for days with multiple saveTags - to avoid overlap, we only take
            %    the savetag with the most trials
            allDates = arrayfun(@(ds) ds.datenum, dc.datasets);
            [uniqueDays, ~, udAssignments] = unique(allDates);

            maskKeep = false(dc.nDatasets, 1);
            nTrials = arrayfun(@(ds) ds.nTrials, dc.datasets);
            for nd = 1:numel(uniqueDays)
                thisDayDatasetInds = find(udAssignments == nd);

                numTrials = nTrials(thisDayDatasetInds);
                [~, whichSetToKeep] = max(numTrials);
                maskKeep(thisDayDatasetInds(whichSetToKeep)) = true;
            end

            dc.filterDatasets(maskKeep);
        end
        
        function filterOutSubsumedSaveTags(dc)
            % remove save tags that are a subset of other save tags
            dates = arrayfun(@(ds) ds.datenum, dc.datasets);
            saveTags = arrayfun(@(ds) ds.saveTagList, dc.datasets, 'UniformOutput', false);
            
            maskKeep = truevec(dc.nDatasets);
            for iD = 1:dc.nDatasets
                maskSubsumed = dates(iD) == dates & cellfun(@(set) all(ismember(saveTags{iD}, set)), saveTags) & maskKeep;
                maskSubsumed(iD) = false; % ignore this dataset
                if any(maskSubsumed)
                    maskKeep(iD) = false;
                end
            end
            
            dc.filterDatasets(maskKeep);
        end
        
        function filterHavingMinimumTrials(dc, minTrials)
            nTrials = cat(1, dc.datasets.nTrials);
            dc.filterDatasets(nTrials >= minTrials);
        end
        
        function filterDatasetsHavingAllConditions(dc)
            if dc.nDatasets == 1
                return;
            end
            
            % unique conditionsacross all datasets 
            allConditions = unique(cat(1, dc.datasets.uniqueConditions));
            
            hasAll = falsevec(dc.nDatasets);
            for iD = 1:dc.nDatasets
                hasAll(iD) = all(ismember(allConditions, dc.datasets(iD).uniqueConditions));
            end
            
            dc.filterDatasets(hasAll);
        end

        function addDataset(dc, ds)
            assert(isa(ds, 'PierreEricLFADS.Dataset'), 'Must be PierreEricLFADS.Dataset instance');
            addDataset@LFADS.DatasetCollection(dc, ds);
        end

        function t = getDatasetInfoTable(dc)
            dc.loadInfo();
            rowNames = arrayfun(@(ds) ds.name, dc.datasets, 'UniformOutput', false);
            date = arrayfun(@(ds) ds.datestr, dc.datasets, 'UniformOutput', false);
            saveTags = arrayfun(@(ds) strjoin(ds.saveTags, ','), dc.datasets, 'UniformOutput', false);
            nChannels = arrayfun(@(ds) ds.nChannels, dc.datasets, 'UniformOutput', true);
            nChannelsHighSNR = arrayfun(@(ds) ds.nChannelsHighSNR, dc.datasets, 'UniformOutput', true);
            nTrials = arrayfun(@(ds) ds.nTrials, dc.datasets, 'UniformOutput', true);

            t = table(date, saveTags, nTrials, nChannels, nChannelsHighSNR, 'RowNames', rowNames);
        end
    end
    
    methods
        function tdSet = loadTrialDataFromDatabase(dc, db, varargin)
            % pierreSettingsEric
            % db = LFADS_PierreExport.Database.loadDb();
            
            p = inputParser();
            p.addOptional('reload', false, @islogical);
            p.addParameter('threshLow', 50, @isscalar);
            p.addParameter('threshHigh', 150, @isscalar);
            p.parse(varargin{:});
            
            if ~isempty(dc.trialDataSet) && ~p.Results.reload
                tdSet = dc.trialDataSet;
                return
            end
            
            if strcmp(dc.context, 'Eric_v5')
                db.loadSource(LFADS_PierreExport.Export.CenterOutExport_v5);
            elseif strcmp(dc.context, 'Eric_v4')
                db.loadSource(LFADS_PierreExport.Export.CenterOutExport_v4);
            elseif strcmp(dc.context, 'Eric_v3')
                db.loadSource(LFADS_PierreExport.Export.CenterOutExport_v3);
            elseif strcmp(dc.context, 'Eric_v2')
                db.loadSource(LFADS_PierreExport.Export.CenterOutExport_v2);
            else
                db.loadSource(PierreDanExport.CenterObutExportv4);
            end
            
            tdSet = cell(dc.nDatasets, 1);
            prog = ProgressBar(dc.nDatasets, 'Loading trialData for each dataset from CenterOutLFADSExport');
            for i = 1:dc.nDatasets
                prog.update(i);
                ds = dc.datasets(i);
                if strcmp(dc.context, 'Eric_v5')
                    st = db.CenterOutLFADSExportv5.match('date', ds.datenum, 'saveTagGroup', ds.saveTags);
                elseif strcmp(dc.context, 'Eric_v4')
                    st = db.CenterOutLFADSExportv4.match('date', ds.datenum, 'saveTagGroup', ds.saveTags);
                elseif strcmp(dc.context, 'Eric_v3')
                    st = db.CenterOutLFADSExportv3.match('date', ds.datenum, 'saveTagGroup', ds.saveTags);
                elseif strcmp(dc.context, 'Eric_v2')
                    st = db.CenterOutLFADSExportv2.match('date', ds.datenum, 'saveTagGroup', ds.saveTags);
                else
                    st = db.CenterOutLFADSExportv4.match('date', ds.datenum, 'saveTag', ds.saveTags);
                end
                td = st.retrieveValue('trialData');
                if ~isempty(td)
                    td = TrialDataConditionAlign(td);
                    
                    if strcmp(dc.context, 'Eric_v2')
                        args = {'smoothing', 31};
                        td = td.dropChannels({'peakSpeed2', 'peakSpeed3', 'handVelocityX', 'handVelocityY', 'handVelocityZ'});
                        td = td.addDifferentiatedAnalogChannel('handX', 'handVelocityX', args{:});
                        td = td.addDifferentiatedAnalogChannel('handY', 'handVelocityY', args{:});
                        td = td.addDifferentiatedAnalogChannel('handZ', 'handVelocityZ', args{:});

                        vx = td.getAnalog('handVelocityX');
                        vy = td.getAnalog('handVelocityY');
                        vz = td.getAnalog('handVelocityZ');

                        s2 = cellfun(@(vx, vy) sqrt(vx.^2 + vy.^2), vx, vy, 'UniformOutput', false);
                        s3 = cellfun(@(vx, vy, vz) sqrt(vx.^2 + vy.^2 + vz.^2), vx, vy, vz, 'UniformOutput', false);

                        % compute speed
                        td = td.addAnalogChannelModifiedInPlace('handVelocityX', 'handSpeed', s2);
                        td = td.addAnalogChannelModifiedInPlace('handVelocityX', 'handSpeed3', s3);

                        % compute peak speed
                        td = td.align('MoveOnsetOnline:+500').truncateAfter('TargetAcquired');
                        % nTrials x nTime x 3 
                        vel = td.getAnalogMultiAsTensor({'handVelocityX', 'handVelocityY', 'handVelocityZ'});
                        s2 = sqrt(sum(vel(:, :, 1:2).^2, 3));
                        s3 = sqrt(sum(vel(:, :, 1:3).^2, 3));
                        max2 = max(s2, [], 2);
                        max3 = max(s3, [], 2);
                        td = td.addScalarParam('peakSpeed2', max2);
                        td = td.addScalarParam('peakSpeed3', max3);
                        td = td.markInvalidTrialsPermanentlyInvalid('Invalid peak speed');
                        
                        % threshold to get better rt at 5% of peak speed in 2D
                        td = td.align('GoCue:GoCue + 600');
                        [speed, speedTime] = td.getAnalog('handSpeed');
                        [speedMat, speedTime] = TrialDataUtilities.Data.embedTimeseriesInMatrix(speed, speedTime);
                        rt = TrialDataUtilities.Data.findThresholdCrossingsLowThenHigh(speedMat', speedTime, max2 * 0.05, max2 * 0.5)';

                        td = td.setParam('rt', rt);
%                         td = td.setEvent('MoveOnsetOnline', rt);
                        td = td.addEvent('Move', rt);

                        % we don't do this in v2 since we hadn't recomputed
                        % the RTs before generating the sequence data
                        
%                         td = td.markTrialsPermanentlyInvalid(isnan(rt), 'Invalid RT');
                        td = td.reset();
                    end
                    
                    tdSet{i} = td;
                else
                    error('Could not locate run %d', i);
                end
            end
            prog.finish();
            
            dc.trialDataSet = tdSet;
        end
        
        function filterDatasets(dc, mask)
            filterDatasets@LFADS.DatasetCollection(dc, mask);
            if ~isempty(dc.trialDataSet)
                dc.trialDataSet = dc.trialDataSet(mask);
            end
        end
    end
end
