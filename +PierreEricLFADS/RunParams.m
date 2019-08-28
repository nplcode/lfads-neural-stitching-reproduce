classdef RunParams < LFADS.RunParams
   properties
       align char = 'MoveOnsetOnline';
       preKeep = 400;
       postKeep = 500; 
       nTrialsKeep = 0;
       
       minDelay = 0;
       maxDelay = Inf;
       
       % for alignment matrix
       pcsKeep uint16 = 8;
       pcTrialAvg logical = true;

       % for running the same model multiple times
       repeatIndex uint16 = 1;
   end
   
   methods
       function par = RunParams()
          % adjust some defaults
          par.spikeBinMs = 10;
          
          par.c_factors_dim = 8;
          par.c_ic_enc_dim = 100;
          par.c_ci_enc_dim = 100;
          par.c_gen_dim = 100;
       end 
        
       function list = getListPropertiesNotAffectingInputDataHash(p)
           list = getListPropertiesNotAffectingInputDataHash@LFADS.RunParams(p);
           list = union(list, {'repeatIndex', 'c_do_readin_bias', 'c_train_readin', 'c_initialize_readout_as'});
       end
   end
end
