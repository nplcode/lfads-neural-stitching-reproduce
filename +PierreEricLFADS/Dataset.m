classdef Dataset < LFADS.Dataset
% A single-day of raw data processed in a common way
    properties(SetAccess=protected)
        nChannelsHighSNR
        uniqueConditions
        
        saveTagList
    end
    
    methods
        function ds = Dataset(collection, relPath)
            ds = ds@LFADS.Dataset(collection, relPath);
        end
        
        function loadInfoFromData(ds, data)
            data = data.data;
            ds.subject = strtok(data.dataset, ' ');
            
            if ismember(ds.collection.context, {'Eric_v4', 'Eric_v5'})
                match = regexp(ds.relPath, 'saveTagGroup_(?<saveTagGroup>[\d,]+)', 'names');
                ds.saveTags = sscanf(match.saveTagGroup, '%d,')';
                match = regexp(ds.relPath, 'saveTag_(?<saveTag>[\d,]+)', 'names');
                ds.saveTagList = sscanf(match.saveTag, '%d,')';
            else
                % find save tag
                match = regexp(data.dataset, 'saveTag(Group)? (?<saveTag>[\d,]+)', 'names');
                ds.saveTags = sscanf(match.saveTag, '%d,')';
            end
            
            % find datestr
            match = regexp(data.dataset, '(?<datestr>\d{4}-\d{2}-\d{2})', 'names');
            ds.datenum  = datenum(match.datestr);
            
            ds.nChannels = data.nUnits;
            ds.nChannelsHighSNR = numel(data.chSNR >= 2);
            ds.nTrials = data.nTrials;
            
            if ismember(ds.collection.context, {'Eric_v2', 'Eric_v3', 'Eric_v4', 'Eric_v5'})
                condField = 'targetDirectionName';
            else
                condField = 'conditionDesc';
            end
            ds.uniqueConditions = unique(data.(condField));
        end
    end
end
