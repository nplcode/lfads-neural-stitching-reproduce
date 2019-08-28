if ~exist('rc', 'var')    
    PierreEricLFADS.Paper.drive_paper_v11_resubmit_fixedReadin;

    db = LFADS_PierreExport.Database.loadDb();
    dc.loadTrialDataFromDatabase(db);

    rc.loadPosteriorMeans();
    rc.loadTrialDataFromDatasetCollection();
    rc.addPosteriorMeansToTrialData();
end

fig_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper_v11_resubmit';
mkdirRecursive(fig_dir);
results_dir = '~/Dropbox/Lab/Projects/LFADS_Stitching/paper_v11_resubmit/results';
mkdirRecursive(results_dir);