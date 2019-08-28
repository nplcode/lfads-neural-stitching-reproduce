% MAKES FIGURE PANEL D

%%%%%%%%%%%%%%%%%%%%
%% Go Cue aligned now
%%%%%%%%%%%%%%%%%%%%%

psetCatGo = PopulationTrajectorySetCrossConditionUtilities.concatenateAlongNewConditionAxis(...
    psetsGo, 'dataset', {ra.datasets.name}', ...
    'equalizeTimeVectors', true, ...
    'aggregateMarks', false, 'aggregateIntervals', false);
%%
psetCatGo = psetCatGo.setConditionAppearanceFn(@PierreEricLFADS.Paper.colorByTarget);

%% Find the CIS projection
[projCIS, psetCatGo_cis] = ProjDPCA_NonOrthogonal.createFromAndProject(psetCatGo, 'lambda', 1e-7, 'nBasesProjPerMarginalization', [1 0 0 0]);

%% Find the first 8 condition dependent modes

psetCatGoTight = psetCatGo.manualSliceOrExpandTimeWindow(0, 480);
dproj = ProjDPCA_NonOrthogonal.createFrom(psetCatGoTight, 'lambda', 1e-7, 'nBasesProjPerMarginalization', [0 8 0 0]);
dproj = dproj.orthonormalize();
psetCatGo_dproj = dproj.projectPopulationTrajectorySet(psetCatGoTight);

%'axesCombineAllMarginalizations', {{{'targetDirectionName'}, {'dataset'}}}

%% Then do jPCA

[jproj, psetCatGo_jproj] = ProjJPCA.createFromAndProject(psetCatGo_dproj, 'normalize', false);
rotProj = StateSpaceProjection.compose(dproj, jproj.filterOutputBases(1:2));

%% Then combine the CIS with the jPCA projection

finalProjGo = StateSpaceProjection.concatenate(projCIS, rotProj);
finalProjGo = finalProjGo.orthonormalize();

[finalCatGo, stats] = finalProjGo.projectPopulationTrajectorySet(psetCatGo);

%%

clf
finalCatGoShow = finalCatGo.manualSliceOrExpandTimeWindow(0, 510);
finalCatGoShow = finalCatGoShow.setConditionAppearanceFn(@PierreEricLFADS.Paper.colorByTarget);
finalCatGoShow = finalCatGoShow.setBasisNames({'CIS', 'jPC1', 'jPC2'});

%%
finalCatGoShow.plotStateSpace('useThreeVector', true, 'zeroShowOnData', true, 'markSize', 15);
axis equal
view(59.2, -58.4)

%% Make the combined plot

figure(2); clf
figSize([9 5]);
g = AutoAxisGrid(1, 2);
g.axisAt(1);

finalCatGoShow.plotStateSpace('alpha', 0.8, 'markOutlineAlpha', 0.5, 'useThreeVector', true, 'zeroShowOnData', true, 'markSize', 12);
axis equal
view(45.6, -52.8)

g.axisAt(2);
finalCatGoShow.plotStateSpace('alpha', 0.8,  'markOutlineAlpha', 0.5, 'useThreeVector', true, 'zeroShowOnData', true, 'markSize', 12);
axis equal
view(90, 0);
g.updatePositions();

%%

saveFigure(fullfile(fig_dir, 'factors.pdf'), 'painters', true, 'upsample', 4);
