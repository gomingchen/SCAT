% Created by Chen Ming, chen_ming@brown.edu on 6/1/2019

% To demonstrate SCAT model searches echoes from lower frequency and only
% analyse echoes with matching lower frequencies with the broadcast
clc;clear;close all;
%% Basic setup
load config_2H_freqHop_telemikecalls_2.mat
load telemike_calls
load call_shizuko_rec2 % the call corresponding to the recording
ts.data = highpass(ts.data, 23E3, ts.fs);
SIML = runBiscatMain2(cfg, ts);
Fs = ts.fs;
nchan = cfg.coch_steps;
nn = 8;

Fc = SIML.coch.Fc;


%% Use xcorr to detect the call after dechirping
[c,lags] = xcorr(ts.data, CALL);
Nn = length(c);
ns = (Nn-1)/2 + 1;
ra = 0.25;
cc = c(ns:end);
lag = lags(ns:end);
figure, plot(lag.*1E3/Fs, cc);
hold on

[yu, ~] = envelope(cc,10,'peak');
plot(lag.*1E3/Fs, yu, 'r', 'LineWidth', 2);

plot([0, 20], ra*max(yu).*[1, 1], 'k', 'LineWidth',2);
thd = max(yu)*ra;
kk = find(yu > thd);

difk = diff(kk);
gk = find(difk>1);
dg = diff(gk);

KK1 = [kk(1); kk(gk+1)];
K = KK1;

%%% Params for half rectify and thresholding process
wavStructParam = struct;
wavStructParam.simStruct = SIML;

wavStructParam.Fs = Fs; % sampling frequency of the broadcast-echoes duo
wavStructParam.callLenForMostFreq = 1.5E-3; % twi, the time when broadcast starts, in seconds;
wavStructParam.callLenForHighFreq = 1.5E-3; % twh
wavStructParam.callLenSpecial = 1.5E-3; % twu
wavStructParam.sepFlag = 1; % if ==1, use fixed Sep, or determine sep with findSep function when all echoes are so close together, i.e., chain echoes
wavStructParam.Sep = 4E-3; % At 4 ms, the broadcast and the echoes are separated
wavStructParam.whenBrStart = 0.3E-3; % br_startSec, the time when broadcast started
wavStructParam.startingThPercent = 5; % th_val, integer, the series of thresholds start from 5 percent.
wavStructParam.th_type = 'const'; % th_type
wavStructParam.SepbwBRand1stEcho = 0.5; % sepearation time between broadcast and 1st echo
wavStructParam.color = 'b';
figure, spectrogram(ts.data,128,120,500,Fs,'yaxis');
figure
threshVec = [4, 5, 6];

for i = 1:3
wavStructParam.NoT = threshVec(i); % which threshold, from 1st to 10th
hold on
[dataStru(i).echoF, firstGap] = linear_separate_window_10thresholds_2H_shizuko(wavStructParam);
end

%%
d = zeros(nchan, length(K)-1, 3);
dbd = zeros(nchan, length(K), 3);
dmo = zeros(nchan.*3, length(K)-1);
dn = zeros(nchan.*3, length(K));
for j = 1:3
[d(:,:,j), dbd(:,:,j)] = categorizeEcho(dataStru(j).echoF, K, Fs, Fc); % dbd, d before dechirp
end

for t = 1:3
dmo(81*(t-1)+1:81*t,1) = d(:,1,t);
dmo(81*(t-1)+1:81*t,2) = d(:,2,t);
dmo(81*(t-1)+1:81*t,3) = d(:,3,t);

dn(81*(t-1)+1:81*t,1) = dbd(:,1,t);
dn(81*(t-1)+1:81*t,2) = dbd(:,2,t);
dn(81*(t-1)+1:81*t,3) = dbd(:,3,t);
dn(81*(t-1)+1:81*t,4) = dbd(:,4,t);
end

dm = zeros(nchan.*3, length(K)-1);
dm(:,1) = dmo(:,1) - min(dmo(:,1)) + 2000;
dm(:,2) = dmo(:,2) - min(dmo(:,2)) + 4000;
dm(:,3) = dmo(:,3) - min(dmo(:,3)) + 6000;
figure
hold on
scatter(zeros(size(Fc)),Fc, 'ob');
for q = 2:3
    if q == 3
        Nbins = 4;
    else
        Nbins = 6;
    end
[a,edges] = histcounts(dm(:,q), Nbins);
scatter(dm(:,q), repmat(Fc', 3,1), 'ob');

    if q == 1
         times = 5;
         midp = 2350;
    elseif q == 2
         times = 20;
         midp = 3950;
    elseif q == 3
         times = 15;
         midp = 5950;
    end
    [st, A] = extendEdges(edges, times, midp);

    for p = 1:length(a)
        patchplot(st+A*(p-1), st+A*p, a(p).*10E3/max(a));
    end

end

%%
figure
scatter(dn(:,1)/Fs.*1E3, repmat(Fc', 3,1), 'ob');
for q = 2:4
    hold on
    scatter(dn(:,q)/Fs.*1E3, repmat(Fc', 3,1), 'ob');
end
    
xlabel('time (ms)');
ylabel('frequency (kHz)');
set(gca, 'YTick', 20E3:10E3:100E3, 'YTickLabel', 20:10:100);
