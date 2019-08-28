function a = colorByTargetHighlightMeanOverDatasets(ci, a, varargin)
    
    targetList = {'Down', 'DownRight', 'Right', 'UpRight', 'Up', 'UpLeft', 'Left', 'DownLeft'};
    centerOutTheta = deg2rad([270 315 0 45 90 135 180 225]);
%     targetCmap = TrialDataUtilities.Color.hslmap(8);
%     evalColorMapAt(TrialDataUtilities.Color.hslmap(1000, 'luminance', 0.4), centerOutTheta, [0 2*pi]);
    
%     targetCmap = [230 25 75; 245 130 48; 255 255 25; 60 180 75; 70 240 240; 0 130 200; 145 30 180; 240 50 230] / 255;
%     targetCMap = TrialDataUtilities.Color.hslmap(8);
    targetCmap = cubehelix(8, 1.09, 2.74, 3.0, 1.34, [0.27 0.80], [0.33 0.74]);

    for iA = 1:ci.nConditions
        cname = ci.names{iA};
%         if strcmp(cname, 'ccm')
%             a(iA).LineWidth = 3;
%             a(iA).Color = [0.2 0.2 0.2];
%             continue;
%         end

        c = ci.conditions(iA);
        a(iA).Color = [0.2 0.2 0.2];

        if ischar(c.targetDirectionName)
            [tf, idx] = ismember(c.targetDirectionName, targetList);
            if tf
                a(iA).Color = targetCmap(idx(1), :);
                a(iA).Alpha = 0.5;
                a(iA).LineWidth = 0.5;
            end
        end
        if strncmp(c.dataset, 'Mean', 4)
            a(iA).LineWidth = 4;
            a(iA).Alpha = 1;
        end
    end
end