function temporal_test
% implements delay-tuned neurons (tapped delay lines) in the spectrogram
% correlation block

clc
clear
close all

% parameters
fs = 1e6;       % "sampling" rate of delay lines
dMax = .003; %0.01;    % delay line length
tMax = 0.1;     % simulation length

% matrix dimensions
M = 80;                % number of CN neurons
N = ceil(dMax*fs);      % number of NLL neurons (per frequency channel)

method = 1;

% use logical sparse matrices for speed
S = sparse([],[],true,M,N);         % initialize a sparse matrix with logical elements
T = 1e3*(0:N-1)/fs;                 % time axis
F = 1e-3./linspace(1/20e3,1/100e3,M)';  % frequency axis

% synthesize pulse and echoes
t0 = 1e3;
t1 = 2e3;
t2 = 2.1e3;

k = 1:M;
S(k,k+t0) = flipud(speye(M))';
S(k,k+t1) = flipud(speye(M))';
S(k,k+t2) = flipud(speye(M))';

% initialize delay line plot as 2D image
fh = figure;
subplot(3,1,1)
raster(T,F,S)
xlabel('Time (ms)')
ylabel('Frequency (kHz)')
title('Delay lines')
colormap(1-gray(2))

%%
% % setup a triggered correlation of the echo with the preceding emission
% % (store the initial volley of spikes per channel)
% for m=1:M
%     D0(m,1) = find(S(m,:),1);
% end

switch method
    case 1
    % method 1 (Saillant, 1993 - Fig. 2a) - uses vertical echo channels along
    % with additional input delays to determine correlation
    
        % find indices of spike locations
        %[I,J] = find(S);
    
        % initialize new "dechirped" sparse representation
        D = sparse([],[],true,M,N);

        % shift each row by delay of first pulse
        for m=1:M
            S0 = find(S(m,:));
            S0 = S0(1) - 1;
            
            % shift row relative to first pulse
            D(m,:) = circshift(S(m,:),[0 -S0]);
            
%             dIdx = sIdx{2} - T0;
%             dIdx(1) = [];           % remove first pulse
%             D(m,dIdx) = true;
        end
    
    case 2
    % method 2 (Saillant, 1993 - Fig. 2b) - uses obliquely oriented echo
    % channel to compensate for timing delay of individual frequency bands
    
        % Not implemented, since we do not assume any linearity between echoes or
        % frequency bands
    
end


%figure;
subplot(3,1,2)
raster(T,F,D)
%xlabel('Time (ms)')
ylabel('Frequency (kHz)')
title('Dechirped delay lines')



%% collapse coincidence detection results across all freqeuncy channels
Dhist = sum(D,1);

%figure;
subplot(3,1,3)
stem(T,Dhist,'.k')
xlabel('Time (ms)')
ylabel('# Impulses')

%% attempted to use xcorr here with decent results, not sure of
% significance of this method though
% x = xcorr(FB.spikes');
% xx=sum(x');
% plot(xx)

%%
% add a dechirping operation by connecting vertical delay lines across
% frequency channels (combines all individual channels into 1 single echo
% channel)


% crude raster plot mode
function raster(X,Y,Z)
nTicks = 8;

[I,J] = find(Z);
plot(X(J),I,'ko',...
            'MarkerFaceColor','k',...
            'MarkerSize',2)

% set Y axis numbering
yInt = floor(length(Y)/nTicks);
yIdx = [(1:yInt:length(Y)-yInt/2) length(Y)];
set(gca,'YTick',yIdx)
set(gca,'YTickLabel',round(Y(yIdx)*10)/10)

axis([X(1) X(end) 1 size(Z,1)])

% image(X,Y,Z); axis xy
% colormap(1-gray(2))

%pcolor(T,F,S)     % use this for rotating planar color image
