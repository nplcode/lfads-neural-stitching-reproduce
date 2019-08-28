function [Tout, model] = kFoldXvalOLE(Tin, K)
% function [Tout, model] = kFoldXvalOLE(Tin, K)


% split into K segments
% use (K-1)/K of the data for train, 1/K of the data for test

allInds = 1:numel(Tin);
for nfold = 0 : K-1
    trainInds = allInds( mod( allInds, K) ~= nfold );
    testInds = setdiff( allInds, trainInds );
    % fit an OLE to the neural data
    model = Analysis.Decoding.fitOLE( Tin(trainInds) );
    % actually run the decoder
    ToutTmp = PierreEricLFADS.Decoding.decodeOLE(model, Tin(testInds));

    % stupid initialization
    if ~exist('Tout','var')
        Tout = ToutTmp;
    end
    Tout(testInds) = ToutTmp;
end
