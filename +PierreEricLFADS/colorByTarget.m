function a = colorByTarget(ci, a, varargin)
    
    targetList = {'Down', 'DownRight', 'Right', 'UpRight', 'Up', 'UpLeft', 'Left', 'DownLeft'};
%     centerOutTheta = deg2rad([270 315 0 45 90 135 180 225]);
%     targetCmap = TrialDataUtilities.Color.evalColorMapAt(TrialDataUtilities.Color.hslmap(1000, 'luminance', 0.4), centerOutTheta, [0 2*pi]);
%     targetCmap = cubehelix(8, 1.09, 2.74, 3.0, 1.34, [0.27 0.80], [0.33 0.74]);
%     targetCmap = circshift([230 25 75; 245 130 48; 182 160 15; 60 180 75; 70 240 240; 0 130 200; 145 30 180; 240 50 230] / 255, -2, 1);
%     targetCmap = circshift([241 86 117; 194 124 38; 132 149 38; 44 166 81; 46 158 146; 47 151 197; 145 115 241; 238 64 213]/255, -1, 1);
    targetCmap = circshift(TrialDataUtilities.Color.hslmap(10, 'luminance', 0.6), -2, 1);
    targetCmap = targetCmap([1 3 5:end], :);
    
    for iA = 1:ci.nConditions
        c = ci.conditions(iA);
        a(iA).Color = [0.2 0.2 0.2];

        if ischar(c.targetDirectionName)
            [tf, idx] = ismember(c.targetDirectionName, targetList);
            if tf
                a(iA).Color = targetCmap(idx(1), :);
                a(iA).Alpha = 1;
                a(iA).LineWidth = 0.5;
            end
        end
%         if strncmp(c.dataset, 'Mean', 4)
%             a(iA).LineWidth = 4;
%             a(iA).Alpha = 1;
%         end
    end
end