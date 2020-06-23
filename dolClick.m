close all; clc; clear;

sandbox

load('dolphin click.mat')
figure, spectrogram(s,128,120,500,Fs,'yaxis');
load config_dolClick_30_150kHz_1kHzStep
%%% Construct wav parameter structure
wavStructParam = struct;
wavStructParam.Fs = Fs; % sampling frequency of the broadcast-echoes duo
wavStructParam.callLenForMostFreq = 2.3E-3; % twi,
wavStructParam.callLenForHighFreq = 1.5E-3; % twh
wavStructParam.callLenSpecial = 2.3E-3; % twu
wavStructParam.sepFlag = 1; % if ==1, use fixed Sep, or determine sep with findSep function when all echoes are so close together, i.e., chain echoes
wavStructParam.whenBrStart = .5E-3; % br_startSec, the time when broadcast started
wavStructParam.startingThPercent = 10; % th_val, integer, the series of thresholds start from 5 percent.
wavStructParam.th_type = 'const'; % th_type
wavStructParam.SepbwBRand1stEchoinSmpls = round(4E-3*Fs); % sepearation time between broadcast and 1st echo, in samples
wavStructParam.color = 'b'; % color setting of the dechirped image.
wavStructParam.ALT = -25; % Coefficient of amplitude latency trading effect, -25 microseconds per dB, must be negative
wavStructParam.NT = 10; % number of thresholds
wavStructParam.ALToff = 1; % flag if ALT is off
wavStructParam.twoHoverlap = [40E3, 55E3]; % the region where two harmonics are overlapped, in kHz
wavStructParam.Necho = 10;
wavStructParam.callLenBr = 2.3E-3; % in seconds
fband = wavStructParam.twoHoverlap;
step = (cfg.coch_fmax - cfg.coch_fmin)/(cfg.coch_steps - 1);
N1 = (fband(1) - cfg.coch_fmin)/step;
N2 = (fband(2) - cfg.coch_fmin)/step;

ts.data = s; %
ts.fs = Fs;
getGlintSpace2(cfg, ts, wavStructParam);


function [Xdelay, delayNN, ip_L, Fc] = getGlintSpace2(cfg, ts, wavPS)
%%% Input the config file and time series of the broadcast-echo duo and get
%%% the delay in microseconds between two glints
figure
SIML = runBiscatMain2_noplot(cfg, ts);
Fc = SIML.coch.Fc;
wavPS.simStruct = SIML;

Necho = 10;  % <-- user input for this example
% 
for t = 1:wavPS.NT
    %%% uncomment this to plot dechirped representations for every
    %%% threshold
    %figure
    %title(sprintf('No. %d threshold',t));
    %%%
    wavPS.NoT = t; % number of thresholds that will be used
    [echoL,~] = linear_separate_window_10thresholds2(wavPS);
    Mat(t,:) = echoL;
    
end

%%% cheat by loading pre-calculated dechirped representations
%load test_dolphinClick
%%%


%%% Find the delay between target and a receiver (could be a ear in
%%% binaural setting or a bat if it's treated as one receiver)

[Ta, d1] = cell2doublematrix(Mat,Necho, wavPS.NT, length(Fc));

ip_L = findnotches_multiechoes_dolphin(Necho, Ta);

%%% Register separate glints if the two echoes from them are apart
% 1. measure width
[index, pos] = measureWidth_dol(ts, 300); % pos here includes the broadcast
% 2. change wave parameter to find the different delays

long_in = setdiff(1:Necho, index);
LA = struct;
comp = pos(9:end)-500-1;
echolength = 500;
for j = 1:length(long_in)
    tj = long_in(j);
    % form a unique broadcast + echo pair
    comp(j) = pos(tj+1) - 500 - 1;
    brd = ts.data(1:pos(1)+ echolength*4);
    echo = ts.data(comp(j)+1: pos(tj+1)+echolength*4);
    ns = struct;
    ns.data = [zeros(500,1); brd; 1.2*echo];
    ns.fs = ts.fs;

    figure
    SIML2 = runBiscatMain2_noplot(cfg, ns);
    wavPS.simStruct = SIML2;
    wavPS.Necho = 2;
    % find the intercepts
    wavPS.startingThPercent = 20;
    wavPS.whenBrStart = 1E-3;
    wavPS.SepbwBRand1stEchoinSmpls = round(6E-3*ns.fs); % sepearation time between broadcast and 1st echo, in samples
    wavPS.callLenForMostFreq = .28E-3; % twi,
    wavPS.callLenForHighFreq = .28E-3; % twh
    wavPS.callLenSpecial = .2E-3; % twu
    
    for t1 = 1:wavPS.NT
        wavPS.NoT = t1;
        [echoL,~] = linear_separate_window_10thresholds_2H_compbio(wavPS);      
        mat(t1,:) = echoL;
    end
    [~, d2] = cell2doublematrix_glint(mat,2, wavPS.NT, length(Fc));
    LA(j).d2 = d2; 
end

%%% cheat by loading pre-calcualted dechirped representations for the
%%% echoes with longer spacing glints -- 8th - 10th echoes
    % load testmat_dolphin_8910
%%%
%% Plot 3D SCAT
A = zeros(length(Fc),Necho);
des = 1:3:length(Fc);
CC = flipud(colormap(winter(length(des))));

for th = 1:8 % threshold 1st to 8th

    for m = 1:Necho
        A(:,m) = d1(:,m,th)/500; 
    end
    Nbin = 5*ones(Necho,1);

figure
Sep = 2;
for w = 1:length(index)

threeDplotfunc_cr_dolphin(A, index(w), des, Nbin, ip_L{w}, Sep, Fc);

end

for q = 1:length(long_in)
    for i = 1:2
       	B = LA(q).d2(:,i,th)/500;
        minD = comp(q)/500+min(LA(q).d2(:,1,th))/500;
        sep = Sep*(q+w-1) - minD;
        add1Hscatter2(comp(q)/500+B, des, sep, Fc, CC, .5);
        delay(i) = addDelayEst_detached_dol(comp(q)/500+B, comp(q)/500+B, sep);
    end
    
    gs = diff(delay);
    str = sprintf('%d \\mus', round(gs*1000));
    text(min(comp(q)/500+B)+sep-.7, Fc(end)/1E3, 20+(q-1)*10, str, 'FontSize', 12, 'FontName', 'arial');
end

adjfig_3dscat_dolphin(th);

end

end