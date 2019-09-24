% Profile_Meddis91

clear
clc
% close all

load 'C:\gaudetteje\src\biscat_signals\synthetic\LFM_tx_only.mat'
% load 'C:\gaudetteje\src\biscat_signals\synthetic\HFM_tx_only.mat'
% load 'C:\gaudetteje\src\biscat_signals\brown_data\JITSIG250kHz.mat'
% load 'C:\gaudetteje\src\biscat_signals\brown_data\JITSIG2MHz.mat'

N = 1;      % # of channels

% Set input level
ts.data = Pascalize(ts.data,80);

% General parameters
model.id = 'Temporal window model (Plack & Oxenham, 1998)';
model.fs = ts.fs; % Sampling frequency
model.ds = 1;     % Outputs downsampled to rate fs/ds

% Model for inner ear
model.cochlea.id='Gammatone filterbank';
model.cochlea.gt.design='gammatone';
model.cochlea.gt.nch=N;
model.cochlea.gt.frange=[44000 55000];

A=AudMod(ts.data,model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Meddis91 parameters for L/M/HSR fibers
%   parameters = [A,B,g,y,l,r,x];
% HSR_params = [20 2500 4000 10.05 1500 6580 66.31];     % [Ca^2+]_thr = 1.4e-11; M = 10; G_Ca_max = 1.6
% MSR_params = [5 300 2000 5.05 2500 6580 66.31];
%MSR_params = [5 300 1000 11.11 1250 16667 250];     % Parameters from Meddis (1988)
MSR_params = [5 300 2000 11.11 1250 16667 250];     % Parameters from Meddis (1988)
% LSR_params = [2 200 500 1.05 3500 6580 66.31];       % [Ca^2+]_thr = 0; M = 10; G_Ca_max = 7.2
aPer = 0.00075;
rPer = 0.0008;
c_r = .55;


% implement meddis91 ANC (output is P[spike]/dt)
% B{1}=ANC_Meddis91(A,ts.fs,HSR_params);
B=ANC_Meddis(A,ts.fs,MSR_params);
% B{3}=ANC_Meddis91(A,ts.fs,LSR_params);

% rand('state',1)
% generate spikes based on P[spike] at each time instant
        Bhat = B;    %tmp
for j=1:300
    Bhat = B(:,1);                % make working copy of P[spike]
    R=rand(size(Bhat));         % init random process
    Chat=zeros(size(Bhat));        % init spike generation
    if aPer
        for ch=1:N
            i=1; idx=1;
            while(i <= length(Bhat) && idx)
                idx = find(Bhat(i:end,ch) > R(i:end,ch), 1);   % iterate over each spike in time, correcting P[spike]
                if isempty(idx), break, end
                Bhat(idx:end,ch) = Ref_Meddis(Bhat(idx:end,1),ts.fs,aPer,rPer,c_r);        % implement meddis91 refractory effects of transmitter available
                Chat(idx,ch) = 1;
            end
        end
        B(:,j)=Bhat;
        C(:,j)=Chat;
    else
        Chat(Bhat > R)=1;
        C(:,j)=Chat;
    end
end

figure;
subplot(3,1,1); plot(ts.time,A); title('Basilar Membrane Movement')
% subplot(4,1,2); plot(B{1}); title('Probability of Spike Occuring (HSR)')
subplot(3,1,2); plot(ts.time,B); title('Probability of Spike Occuring (MSR)')
% subplot(4,1,4); plot(B{3}); title('Probability of Spike Occuring (LSR)')

% figure;
% subplot(3,1,1); plot(C{1}); set(gca,'YLim',[-.1 1.1]); title('Auditory Nerve Output (HSR)')
subplot(3,1,3); plot(ts.time,C); set(gca,'YLim',[-.1 1.1]); title('Auditory Nerve Output (MSR)')
% subplot(3,1,3); plot(C{3}); set(gca,'YLim',[-.1 1.1]); title('Auditory Nerve Output (LSR)')

%BMMplot(A,[1:4],ts.fs,[1:4]);