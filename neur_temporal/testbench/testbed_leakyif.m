%% testbed for leaky integrate and fire neurons with synaptic inputs

clc
clear
close all

rand('state',3)     % for debugging

% simulation specific parameters
tLen = 1;
fs = 1e4;
t = linspace(0,tLen,ceil(tLen*fs))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate synaptic input spike trains
% spike generator parameters
N = 100;                        % Number of spike train processes
M = 4;                          % Number of postsynaptic neurons
r0 = 150;                        % Poisson process mean spiking rate (Hz)
tau_ref = 0.001;                % Exponential refractory time constants (sec)

% Assign synaptic weights to postsynaptic neurons
W = .5*ones(M,N);                  % Input synaptic weights
% wIdx = [0 2 ; 1 3 ; 2 4 ; 3 5] * 100/5 + [ones(M,1) zeros(M,1)]
W(1,1:40) = 0;
W(2,21:60) = 0;
W(3,41:80) = 0;
W(4,61:100) = 0;
% W = rand(M,N);                  % Randomize synaptic weights
R = zeros(M);                   % Absent recurrent synaptic weights
% R = 1 - eye(M);                     % Strong excitatory recurrent synaptic weights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate N spike trains for refractory period
rcos = 100*(1+cos(2*pi*10*t));
tic
for n=1:N
    % Poisson generator with constant firing rate and refractory period
    S(:,n) = genspikes(r0,tau_ref,length(t),fs,'bernoulli');
    
    % Poisson generator with sinusoidally varying firing rate
%     S(:,n) = genspikes(rcos,0,length(t),fs,'bernoulli');

end
toc


%% plot spike trains
if N <= 50
    S_hat = S + ones(length(S),1) * (0:2:2*N-1);
    figure;
    plot(t,S_hat,'k')
    axis([0 tLen -0.1 2*N+0.1])
    grid on;
end

%% Plot ISI histogram of synaptic input
Iin = genisi(S(:),fs);        % calculate interspike-interval firing rate
Iin = Iin(:);
CvIn = std(Iin)/mean(Iin);
figure; histpdf(Iin,100,[0 max(Iin)]);
title('ISI Histogram for all synaptic input spike trains')
a = axis;
text(.6*a(2), .75*a(4), sprintf('C_v = %.2f',CvIn));
clear a


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate IF neurons

% internal model parameters (assume homogeneous neurons) - copied from
% Song, Miller, and Abbott
NN.Vth = -54e-3;                    % neural spiking threshold
NN.Vrest = -70e-3;                  % rest potential
NN.Vreset = -80e-3;                 % reset potential
NN.tau_m = 10e-3;                   % time constant of membrane potential (leakiness factor)
NN.tau_ex = 5e-3;                   % time constant of excitatory synapse
NN.tau_in = 5e-3;                   % time constant of inhibitory synapse
NN.tau_re = 5e-3;                   % time constant of recurrent synapse

% generate spikes using the leaky integrate and fire model
tic
[V,g_ex,g_in,g_re] = leakyifsyn(t, NN, S, W, R, fs);
toc

%% calculate statistics on V and ISI
V_mu = mean(V,1);
V_std = std(V,1);
Iout = genisi(V,fs);
Cv_out = std(Iout)/mean(Iout);

%% Plot neuron membrane potential along with external stimulus vs. time
figure;
plot(t,V)
grid on; hold on;
plot([0 tLen],[NN.Vth NN.Vth],'--r')
title(sprintf('Integrate and Fire Neuron Membrane Potential (N=%d, M=%d)', N, M))
xlabel('Time (sec)')
ylabel('Potential (V)')
set(gca,'YLim',[NN.Vreset*1.2 0.12])

%% Plot neuron input transconducances (g_ex, g_in)
figure;
plot(t,g_ex,t,g_in)
grid on;
title(sprintf('Excitatory and Inhibitory Transconductances (N=%d, M=%d)', N, M))
xlabel('Time (sec)')
ylabel('Transconductance (dimensionless)')

%% Plot neuron recurrent transconducances (g_re)
figure;
plot(t,g_re)
grid on;
title(sprintf('Recurrent Excitatory Transconductances (N=%d, M=%d)', N, M))
xlabel('Time (sec)')
ylabel('Transconductance (dimensionless)')

%% Plot ISI histogram of output neurons
if ~isempty(Iout)
    figure;
    histpdf(Iout,100,[0 max(Iout)])
    title('ISI Histogram of CN Spikes')
    a = axis;
    text(.6*a(2), .75*a(4), sprintf('C_v = %.2f',CvIn));
    clear a
end

%% Plot autocorrelation of output neurons
Sout = zeros(size(V,1),1); Sout(V(:,1)>0) = 1;
if sum(Sout) > 3
    [A,M] = autocorrspike(Sout(:,1),fs,0.2);
    figure;
    bar(M*1e3,A,1);
    title('Autocovariance of Bushy cell response')
    xlabel('Time (ms)')
    grid on
end
