function cn = runCochNuc(neur, ts, Fc, cfg)
% RUNCOCHNUC  Generates spikes in a Bushy cell feedforward/recurrent network

% Author:  Jason Gaudette
% Date:    5/08/2009
% Version: 0.3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feedforward/Recurrent Network of Bushy Cells: Leaky Integrate and Fire Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign user parameters
M = numeric(cfg.cn_ntot);
Nsyn = numeric(cfg.cn_nsyn);
Ntot = size(neur.spikes,2);
wRec = numeric(cfg.cn_wrec);
nRec = numeric(cfg.cn_nrec);

% Assign synaptic weights to postsynaptic neurons
% W = rand(M,Ntot);                  % Randomize synaptic weights
W = zeros(M,Ntot);
Ninc = (Ntot-Nsyn)/(M-1);
idx = 0;
assert(Ntot >= Nsyn, 'Number of synapses per cell (%d) must be greater than total number of channels (%d)',Nsyn,Ntot)
for m=1:M
    idx = (1:Nsyn) + round(Ninc*(m-1));
    W(m,idx) = 1;
end

% Assign recurrent synaptic weights
R = wRec * (1 - eye(M));                 % Strong excitatory recurrent synaptic weights
%r = wRec * (1 - eye(nRec));        % apply only to nearest N neighbors


%% additional LIF parameters - need to add this to leakyifsyn.m
% g_exDel = 0.015;            % define incremental value of g_ex after each presynaptic AP
% g_inDel = 0.05;             % define incremental value of g_in after each presynaptic AP
% g_reDel = 0.05;              % define incremental value of g_re after each recurrent AP
% E_ex = 0;                   % excitatory reversal potential of 0mV
% E_in = -0.07;               % inhibitory reversal potential of -70mV



%% set parameters for homogenous LIF neurons

% internal model parameters (assume homogeneous neurons)
cfg.cn_lifparams = cfg.cn_lifparams * 1e-3;
cn.NN.Vth = cfg.cn_lifparams(1);        % neural spiking threshold
cn.NN.Vrest = cfg.cn_lifparams(2);      % rest potential
cn.NN.Vreset = cfg.cn_lifparams(3);     % reset potential
cn.NN.tau_m = cfg.cn_lifparams(4);      % time constant of membrane potential (leakiness factor)
cn.NN.tau_ex = cfg.cn_lifparams(5);     % time constant of excitatory synapse
cn.NN.tau_in = cfg.cn_lifparams(6);     % time constant of inhibitory synapse
cn.NN.tau_re = cfg.cn_lifparams(7);     % time constant of recurrent synapse
cn.NN.g_exDel = cfg.cn_lifparams(8);    % define incremental value of g_ex after each EPSP
cn.NN.g_inDel = cfg.cn_lifparams(9);    % define incremental value of g_in after each IPSP
cn.NN.g_reDel = cfg.cn_lifparams(10);   % define incremental value of g_re after each recurrent PSP

if 0        % for test purposes only
    cn.NN.Vth = -54e-3;
    cn.NN.Vrest = -70e-3;
    cn.NN.Vreset = -80e-3;
    cn.NN.tau_m = 10e-3;
    cn.NN.tau_ex = 5e-3;
    cn.NN.tau_in = 5e-3;
    cn.NN.tau_re = 5e-3;
    cn.NN.g_exDel = 15e-3;
    cn.NN.g_inDel = 5e-3;
    cn.NN.g_reDel = 5e-3;
end


%% simulate membrane potentials using the leaky integrate and fire model
[cn.V, cn.g_ex, cn.g_in, cn.g_re] = leakyifsyn(ts.time, cn.NN, neur.spikes, W, R, ts.fs);

%% generate discrete spike trains using resulting potential
cn.spikes = zeros(size(cn.V));
cn.spikes(cn.V > 0) = 1;
cn.time = ts.time;
cn.fs = ts.fs;

%% calculate statistics on V and ISI
cn.V_mu = mean(cn.spikes,1);
cn.V_std = std(cn.spikes,1);
cn.ISI = genisi(cn.spikes,ts.fs);
cn.Cv = std(cn.ISI)/mean(cn.ISI);

function param = numeric(param)
% converts string input parameters to numeric type
if ischar(param)
    if numel(param) > 1
        param = str2num(param);
    else
        param = str2double(param);
    end
end
