function S = leakyifsyn_testbed
%% testbed for integrate and fire network

% setup parameters
fs = 1e5;                       % Sampling rate (Hz)
tLen = 1;                       % Epoch duration (sec)

% spike parameters
N = 20;                          % Number of spike train processes
r0 = 100;                       % Poisson process mean spiking rate (Hz)
tau_ref = .01;                  % Exponential refractory time constants (sec)

% histogram parameters
bWidth = 5e-3;                  % ISI histogram bin width (sec)
bRange = [0 0.1];               % ISI histogram data range (sec)

%%
rand('state',1)

% Generate N spike trains for refractory period of 10ms
for n=1:N
    S(:,n) = genspikes(r0,0.01,tLen,fs,'bernoulli');
    I{n} = genisi(S(:,n),fs);
    CoefVarISI(n) = std(I{n})/mean(I{n});     % calculate coefficient of variation
end

% % plot histogram of ISI
% figure;
% histpdf(I{1},ceil(diff(bRange)/bWidth),bRange);
% grid on;
% 
% % plot spike train using raster plot
% 
% Shat = S + ones(length(S),1) * (0:2:2*N-1);
% T = linspace(0,tLen,length(S))';
% figure;
% plot(T,Shat)
% axis([0 tLen -0.1 2*N+0.1])
% grid on;
% 
% 
% figure;
% [Sx,Sy]=find(S);
% plot(Sx,Sy,'ok')
