close all; clc; clear;

sandbox
[s, Fs] = audioread('new 2 glint LP series.wav');


figure, spectrogram(s,128,120,500,Fs,'yaxis');
load config_halfkHzStep
%%% Construct wav parameter structure
wavStructParam = struct;
wavStructParam.Fs = Fs; % sampling frequency of the broadcast-echoes duo
wavStructParam.callLenForMostFreq = 1.5E-3; % twi,
wavStructParam.callLenForHighFreq = 1.5E-3; % twh
wavStructParam.callLenSpecial = 1.5E-3; % twu
wavStructParam.sepFlag = 1; % if ==1, use fixed Sep, or determine sep with findSep function when all echoes are so close together, i.e., chain echoes
wavStructParam.whenBrStart = .5E-3; % br_startSec, the time when broadcast started
wavStructParam.startingThPercent = 3; % th_val, integer, the series of thresholds start from 5 percent.
wavStructParam.th_type = 'const'; % th_type
wavStructParam.SepbwBRand1stEchoinSmpls = round(5E-3*Fs); % sepearation time between broadcast and 1st echo, in samples
wavStructParam.color = 'b'; % color setting of the dechirped image.
wavStructParam.ALT = -25; % Coefficient of amplitude latency trading effect, -25 microseconds per dB, must be negative
wavStructParam.NT = 10; % number of thresholds
ts.data = s(1:24E3); %
ts.fs = Fs;
getGlintSpace2_cr(cfg, ts, wavStructParam);



function getGlintSpace2_cr(cfg, ts, wavPS)
%%% Input the config file and time series of the broadcast-echo duo and get
%%% the delay in microseconds between two glints
figure
SIML = runBiscatMain2_noplot(cfg, ts);
Fc = SIML.coch.Fc;
wavPS.simStruct = SIML;

Necho = 7;

%%% comment this for loop and load test_new_CR_sequence
for t = 1:wavPS.NT

    wavPS.NoT = t; % number of thresholds that will be used
    [echoL,~] = linear_separate_window_10thresholds(wavPS);
    %d1(t,:,) = firstGapL;
    Mat(t,:) = echoL;
    title(sprintf('No. %d threshold',t));
    
end

%%% Cheat by uncomment the line below to skip the calculation above
    %load test_new_CR_sequence


%%% Find the delay between target and a receiver (could be a ear in
%%% binaural setting or a bat if it's treated as one receiver)

[Ta, d1] = cell2doublematrix(Mat,Necho, wavPS.NT, length(Fc));

ip_L = findnotches_multiechoes(Necho, Ta);

ip_ph = findphantomNotches(Necho, Ta, ip_L, 161, 10); % find phantom "notches" that are activated because of ALT


%% Plot 3D SCAT
Fc = 20E3:2E3:100E3;
des = 1:3:161;
opind = 4*ones(Necho,1);
A = zeros(161,Necho);
for th = 2:8

    for m = 1:Necho
        A(:,m) = d1(:,m,th)/500; 
    end
    Nbin = 5*ones(Necho,1);



figure
intvl = .5E3;
intvl_thin = 2E3;
Sep = 2;
Nn_t = zeros(Necho, 40);
for w = 1:Necho
ph_null = ip_ph(th).ft;
LAT = ip_ph(th).lat;
Nn_t(w,:) = threeDplotfunc_cr(A, w, des, intvl,intvl_thin, Nbin, ip_L{w}, Sep, ph_null(:,w), LAT);

end

adjfig_3dscatCR(th);

end

end
