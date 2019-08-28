%% Load data

fig_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper_v11_resubmit/';
results_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper_v11_resubmit/results';
export_file = fullfile(results_dir, 'exported_single_trial_rosette_for_video.mat');
export = load(export_file);

% flip dim 1 to make CIS increasing
% flds = fieldnames(export);
% for f = 1:numel(flds)
%     for i = 1:numel(export.(flds{f}))
%         p = export.(flds{f}){i};
%         p.dataByTrial{1} = -p.dataByTrial{1};
%         p.dataMean{1}(1, :, :) = -p.dataMean{1}(1, :, :);
%         export.(flds{f}){i} = p;
%     end
% end

%%


%% Setup the figure and video

figh = figure(1); clf;
% figh.ToolBar = 'none';
% figh.MenuBar = 'none';
figSize([18 3]);

az = 59.4;
el = -58.4;
C = 7;
g = AutoAxisGrid(1, C);
for i = 1:C
    g.axisAt(i);
  
    p = export.factors_stitch_indiv_rosette_tight{i};
    p.dataByTrial{1} = -p.dataByTrial{1};
    p.dataMean{1}(1, :, :) = -p.dataMean{1}(1, :, :);
    plotSingleTrials(p, 0, 510);
    view(az, el);
%     if i == 1
%         xlim([-0.5245    0.2755]);
%         ylim([-0.7200    0.5800]);
%         zlim([-0.6052    0.6948]);
%     else
%         xlim([-0.4 0.4]);
%         ylim([-0.5 0.8]);
%         zlim([-0.5 0.8]);
%     end

end

g.spacing_y = [0.2 -0.2];
g.spacing_x = [0 -2 -0.5];
g.updatePositions();

axh = axes('Units', 'normalized', 'Color', 'none', 'HitTest', 'off');
axh.OuterPosition = [0 0 1 1];
xlabel(axh, 'CIS');
ylabel(axh, 'jPC1');
zlabel(axh, 'jPC2');
axis(axh, 'off');
view(axh, az, el);
axis(axh, 'equal');

axhText = axes('Units', 'normalized', 'Color', 'none', 'HitTest', 'off');
axhText.Position = [0 0 1 1];
xlim(axhText, [0 1]);
ylim(axhText, [0 1]);
axis(axhText, 'off');
for c = 1:C
    pos = g.handles{1, c}.Position;
    h = text(pos(1) + pos(3)/2, 0.9, sprintf('Session %d', c), 'Parent', axhText, 'Background', 'none', 'HorizontalAlignment', 'center');
end

tv = ThreeVector(axh, false);
tv.background = false;
tv.vectorLength = 0.5;
tv.textVectorNormalizedPosition = 2;
tv.lineColor = [0.4 0.4 0.4];
tv.fontColor = [0.1 0.1 0.1];
tv.update();

% adjust jPC2 label manually
tv.ht(3).Position = [0.0160    0.3286         0]; 

%%

% saveFigure(fullfile(fig_dir, 'singleTrialRosettes.pdf'), 'painters', false, 'upsample', 10, 'resolution', 720);

saveFigure(fullfile(fig_dir, 'singleTrialRosettes.pdf'), 'painters', true);

%%

function plotSingleTrials(p, tMin, tMax)
    tMask = p.tvecDataByTrial{1} >= tMin & p.tvecDataByTrial{1} <= tMax;
    
    % plot go cues
    sx = p.dataByTrial{1}(:, find(tMask, 1));
    sy = p.dataByTrial{2}(:, find(tMask, 1));
    sz = p.dataByTrial{3}(:, find(tMask, 1));
    
    scatter3(sx, sy, sz, 10, 'Clipping', 'off', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerEdgeAlpha', 0.5, 'LineWidth', 1);
    hold('on');
    
    % plot trajectories
    for c = 1:p.nConditions
        sx = p.dataByTrial{1}(p.trialLists{1, c}, tMask)';
        sy = p.dataByTrial{2}(p.trialLists{1, c}, tMask)';
        sz = p.dataByTrial{3}(p.trialLists{1, c}, tMask)';
        
        color = p.conditionDescriptor.appearances(c).Color;
        
        plot3(sx, sy, sz, '-', 'LineWidth', 0.5, 'Color', [color 0.3], 'Clipping', 'off');box 
    end
    hold('off');
    xlabel('CIS'); ylabel('jPC1'); zlabel('jPC2');
    
    axis equal
    axis tight
    axis manual
    axis off;
    
%     tv = ThreeVector(gca, true);
% %     tv.flipAxis(1) = true;
% %     tv.flipAxis(2) = true;
%     tv.fontSize = 8;
%     tv.vectorLength = 1.5;
%     tv.update();
end
