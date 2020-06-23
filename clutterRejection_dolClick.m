close all; clc; clear;
% 
sandbox
load('dolphin click clutter rej.mat');

figure, spectrogram(s,128,120,500,Fs,'yaxis');
load config_dolClick_30_150kHz_1kHzStep
%%% Construct wav parameter structure
wavStructParam = struct;
wavStructParam.Fs = Fs; % sampling frequency of the broadcast-echoes duo
wavStructParam.callLenForMostFreq = 1E-3; % twi,
wavStructParam.callLenForHighFreq = 1E-3; % twh
wavStructParam.callLenSpecial = 1E-3; % twu
wavStructParam.sepFlag = 1; % if ==1, use fixed Sep, or determine sep with findSep function when all echoes are so close together, i.e., chain echoes
wavStructParam.whenBrStart = .5E-3; % br_startSec, the time when broadcast started
wavStructParam.startingThPercent = 5; % th_val, integer, the series of thresholds start from 5 percent.
wavStructParam.th_type = 'const'; % th_type
wavStructParam.SepbwBRand1stEchoinSmpls = round(2E-3*Fs); % sepearation time between broadcast and 1st echo, in samples
wavStructParam.color = 'b'; % color setting of the dechirped image.
wavStructParam.ALT = -25; % Coefficient of amplitude latency trading effect, -25 microseconds per dB, must be negative
wavStructParam.NT = 10; % number of thresholds

ts.data = s; %
ts.fs = Fs;
getGlintSpace2_cr2(cfg, ts, wavStructParam);



function getGlintSpace2_cr2(cfg, ts, wavPS)
%%% Input the config file and time series of the broadcast-echo duo and get
%%% the delay in microseconds between two glints
figure
SIML = runBiscatMain2_noplot(cfg, ts);
Fc = SIML.coch.Fc;
wavPS.simStruct = SIML;

Necho = 2; % <-- user input

%%% comment this for loop and load test_Dolphin_Clutter_Rej.mat
for t = 1:wavPS.NT
    %figure
    wavPS.NoT = t; % number of thresholds that will be used
    [echoL,~] = linear_separate_window_10thresholds2(wavPS);
    %d1(t,:,) = firstGapL;
    Mat(t,:) = echoL;
    title(sprintf('No. %d threshold',t));
    
end

%%% Cheat by uncomment the line below to skip the calculation above
    %load test_Dolphin_Clutter_Rej


%%% Find the delay between target and a receiver (could be a ear in
%%% binaural setting or a bat if it's treated as one receiver)

[Ta, d1] = cell2doublematrix2(Mat,Necho, wavPS.NT, length(Fc),500);
ip_L = findnotches_multiechoes(Necho, Ta);
ip_ph = findphantomNotches2(Necho, Ta, Fc, ip_L, cfg.coch_steps, 10); % find phantom "notches" that are activated because of ALT


%% Plot 3D SCAT

des = 1:3:length(Fc);
A = zeros(cfg.coch_steps,Necho);
for th = 1:9

    for m = 1:Necho
        A(:,m) = d1(:,m,th)/500; 
    end
    Nbin = 5*ones(Necho,1);

figure
intvl = 1E3;
intvl_thin = 2E3;
Sep = 2;
%Nn_t = zeros(Necho, 40);
for w = 1:Necho
ph_null = ip_ph(th).ft;
LAT = ip_ph(th).lat;
threeDplotfunc_cr_dol(A, Fc, w, des, intvl,intvl_thin, Nbin, ip_L{w}, Sep, ph_null(:,w), LAT);

end

adjfig_3dscat_dolCR(th);
end

end
