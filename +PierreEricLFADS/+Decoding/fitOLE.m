function model = fitOLE( data , opts )
% function model = fitOLE( data, opts )
%
% CP: script to fit simple Optimal Linear Estimator on neural and
%   kinematic data. Based on code from Jon Kao.
%
% data to be a trial-ized struct
%   data(n).T = number of time bins
%   data(n).X(3:4,1:T) should be the X and Y velocities
%   data(n).X(5,1:T) should be 1 x T
%   data(n).Z = neural data, dims of numChannels x T
%   data(1).dt = bin time in ms
%
% if opts.channels is set, obs model will only be fit on those channels

opts.foo = false;
if ~isfield(opts,'channels')
    opts.channels = 1:size(data(1).Z,1);
end
channels = opts.channels;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Fit Optimal Linear Estimator
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % djoshea adding for trial masking
  if isfield(data, 'valid')
    data = data([data.valid]);
  end
  
  % add a row of 1's in the neural data.  The kinematic data
  % should only consider the 2d position and velocity.

  idx2D = [1 2 3 4];
  Z     = [data.Z];

  % append a row of 1's to Z.
  Z     = [Z(channels,:); ones(1, size(Z,2))];
  X     = [data.X];
  X     = X(idx2D, :);

  % solve the Wiener problem
  Lw    = mrdivide(X, Z);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Assign model parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%

  model.Lw   = Lw;
  model.dt  = data(1).dt;
