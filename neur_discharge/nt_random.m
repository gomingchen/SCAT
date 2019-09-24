function neur = nt_random(ts, cfg)
% NT_RANDOM  Generates spikes based on a random Poisson process and predefined
% (time-dependent) probability of spiking
%
% neur = NT_RANDOM(TS,CFG)

% define user parameters
rPer = numeric(cfg.neur_rand_rper)*1e-3;
aMod = numeric(cfg.neur_rand_amod);
Ntot = cfg.coch_steps * numeric(cfg.neur_rand_fanout);
T = length(ts.time);

% set spike generation firing rate per channel
switch cfg.neur_rand_mode
    case 1              % constant rate
        % Poisson generator with constant firing rate
        fRate = numeric(cfg.neur_rand_rate);
    case 2              % sinusoidal rate
        % Poisson generator with sinusoidally varying firing rate
        fRate = numeric(cfg.neur_rand_rate)*(1+cos(2*pi*aMod*ts.time));
    case 3              % bypass
        fRate = 0;
    otherwise
        fRate = 0;
        warning('NeurTran:MODE','Unknown random spike generation mode, letting fRate to 0')
end

% Generate random poisson spike trains for each channel
if (rPer > 0)
    neur.spikes = genspikes(fRate,rPer,[T Ntot],ts.fs,'bernoulli');
else
    neur.spikes = genspikes(fRate,0,[T Ntot],ts.fs,'bernoulli');
end

% Return firing rate as probability of spike
if length(fRate) == 1
    neur.prob = fRate*ones(T,Ntot);
else
    neur.prob = fRate*ones(1,Ntot);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param = numeric(param)
% converts string input parameters to numeric type

if ischar(param)
    if numel(param) > 1
        param = str2num(param); %#ok<ST2NM>
    else
        param = str2double(param);
    end
end
