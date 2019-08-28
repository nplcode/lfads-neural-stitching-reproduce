function Tout = decodeOLE(M, T)
% function out = decodeOLE(M, T)
% 
% CP: script to decode kinematics from neural data using simple Optimal Linear Estimator. Based on code from Jon Kao.
%
% M is model
% T is Tstruct

% alpha is a parameter for mixing in position decoding
%   in Jon's paper this was 0.975. Here we will use 1

alpha = 1;

for ntr = 1:length(T)
    if isfield(T, 'valid')
        Tout(ntr).valid = T(ntr).valid; %#ok<*AGROW>
    end
    Tout(ntr).xk = nan( size(T(ntr).X ));
    % give the starting kinematics for every trial
    Tout(ntr).xk(1:2, 1) = T(ntr).X( 1:2, 1);
    neural = T(ntr).Z(:, 1);
    neural = [neural; 1];
    % decode
    z  = M.Lw * neural;
    Tout(ntr).xk(3:4, 1) = z(3:4);
    Tout(ntr).xk(5, 1) = 1;

    for nt = 2 : size( T(ntr).X, 2)
        neural = T(ntr).Z(:, nt);
        neural = [neural; 1];
        z  = M.Lw * neural;

        Tout(ntr).xk(1:2, nt) = (1-alpha)*z(1:2) + ...
            alpha * ( Tout(ntr).xk(1:2, nt-1) + z(3:4)*M.dt / 1000); % MODIFIED TO USE VELOCITIES PER SEC
        Tout(ntr).xk(3:4, nt) = z(3:4);
        Tout(ntr).xk(5, nt)   = 1;
    end
    
    Tout(ntr).T = T(ntr).T;
    Tout(ntr).dt = T(ntr).dt;
    Tout(ntr).time = T(ntr).time;
    Tout(ntr).condition = T(ntr).condition;
end
