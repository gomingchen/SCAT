% This program is to demonstrate the amplitude-latency trading effect of
% auditory-system-inspired SCAT model

clc;clear;close all;
%%% Load configuration file
load config_halfkHzStep.mat

%%% Generate pulse-echo pair with two glints
dist2tg = 1; % distance between sound source and target, in meters
glintDelay = 30; % the delay between two glints in one target, in microsecond
ts = generate_sigs_with_delay_multiglints(dist2tg, glintDelay);
Fs = ts.fs;

%%% Construct wav parameter structure
wavStructParam = struct;
wavStructParam.Fs = Fs; % sampling frequency of the broadcast-echoes duo
wavStructParam.callLenForMostFreq = 1.5E-3; % twi,
wavStructParam.callLenForHighFreq = 1.5E-3; % twh
wavStructParam.callLenSpecial = 1.5E-3; % twu
wavStructParam.sepFlag = 1; % if ==1, use fixed Sep, or determine sep with findSep function when all echoes are so close together, i.e., chain echoes
wavStructParam.Sep = 8E-3; % At 8 ms, the broadcast and the echoes are separated
wavStructParam.whenBrStart = 0.3E-3; % br_startSec, the time when broadcast started
wavStructParam.startingThPercent = 5; % th_val, integer, the series of thresholds start from 5 percent.
wavStructParam.th_type = 'const'; % th_type
wavStructParam.SepbwBRand1stEchoinSmpls = 8*500; % sepearation time between broadcast and 1st echo, in samples
wavStructParam.color = 'r'; % color setting of the dechirped image.
wavStructParam.ALT = -25; % Coefficient of amplitude latency trading effect, -25 microseconds per dB, must be negative
%%% Find the dechirped image
[deL1, ipL, Fc] = getGlintSpace2(cfg, ts, wavStructParam); % in microseconds

%%% Plot the notches onto the network
hearcellspread_empty_toomanychannels(Fc);
u = 0;
if length(ipL)>1
    [fL, sL] = hearcellspreaduniversal_2e_nulls_half(Fc, ipL, 'r','-r');
    u = 1;
end
title('Triangular Network');


