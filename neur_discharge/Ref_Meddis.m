function P = Ref_Meddis(P, fs, varargin)
% Ref_Meddis  Implementation of Meddis' (1991) neural refractory period
%
% P = Ref_Meddis(P, Fs, rPer)
%

% define default model parameters
tau_abs = 0.75e-3;
tau_rel = 0.8e-3;
c_r = 0.55;

% parse optional user paramaters
switch nargin
    case 2
    case 3
        tau_abs = varargin{1};
    case 4
        tau_abs = varargin{1};
        tau_rel = varargin{2};
    case 5
        tau_abs = varargin{1};
        tau_rel = varargin{2};
        c_r = varargin{3};
    otherwise
        error('Incorrect number of input parameters used')
end

% set absolute refractory period (if non-zero)
if tau_abs
    idx = min(length(P),round(fs*tau_abs));
    P(1:idx) = 0;
else
    idx = 1;
end

% set relative refractory period
P(idx:end) = P(idx:end) .* (1 - c_r * exp(-(0:length(P)-idx)'/(tau_rel*fs)));