close all; clc; clear;

sandbox
[s, Fs] = audioread('equalized bat FM signals with Uday 2 glint series.wav');


load config_halfkHzStep
cfg.coch_bw = 2000;

%%% Construct wave parameter structure
wavStructParam = struct;
wavStructParam.Fs = Fs; % sampling frequency of the broadcast-echoes duo
wavStructParam.callLenForMostFreq = 2.3E-3; % twi,
wavStructParam.callLenForHighFreq = 1.5E-3; % twh
wavStructParam.callLenSpecial = 2.3E-3; % twu
wavStructParam.callLenBr = 2.3E-3; % call length of broadcast
wavStructParam.sepFlag = 1; % if ==1, use fixed Sep, or determine sep with findSep function when all echoes are so close together, i.e., chain echoes
wavStructParam.whenBrStart = .5E-3; % br_startSec, the time when broadcast started
wavStructParam.startingThPercent = 20; % th_val, integer, the series of thresholds start from 5 percent.
wavStructParam.th_type = 'const'; % th_type
wavStructParam.SepbwBRand1stEchoinSmpls = round(4E-3*Fs); % sepearation time between broadcast and 1st echo, in samples
wavStructParam.color = 'b'; % color setting of the dechirped image.
wavStructParam.ALT = -25; % Coefficient of amplitude latency trading effect, -25 microseconds per dB, must be negative
wavStructParam.NT = 10; % number of thresholds
wavStructParam.ALToff = 1; % flag if ALT is off
wavStructParam.twoHoverlap = [40E3, 55E3]; % the region where two harmonics are overlapped, in kHz
wavStructParam.Necho = 10;
fband = wavStructParam.twoHoverlap;
step = (cfg.coch_fmax - cfg.coch_fmin)/(cfg.coch_steps - 1);
N1 = (fband(1) - cfg.coch_fmin)/step;
N2 = (fband(2) - cfg.coch_fmin)/step;

wavStructParam.twoHarmonicNNs = [N1 N2]; % the indices for the overlap of the two harmonics
ts.data = s;
ts.fs = Fs;
getGlintSpace2H(cfg, ts, wavStructParam);



function getGlintSpace2H(cfg, ts, wavPS)
%%% Input the config file and time series of the broadcast-echo duo and get
%%% the delay in microseconds between two glints
figure
SIML = runBiscatMain2_noplot(cfg, ts);
Fc = SIML.coch.Fc;
wavPS.simStruct = SIML;
Necho = wavPS.Necho;

matEcho = {};
matB = {};
Nchan = length(Fc);
de = struct; % the dechirped 1st & 2nd harmonics
echolength = 500; %  number of samples

% Find out the echoes whose glint spacing is less than 500 us to be
% analyzed by SCAT

plotflag = 0; % 0: no plots for intercepts and dechirped representations

g = 0;
for t = 1:wavPS.NT
    %%% uncomment to generate dechirped 
        %figure
        %title(sprintf('No. %d threshold',t));
    %%% 
    g = g + 1;
    hold on
    wavPS.NoT = t; % number of thresholds that will be used
    [echoL,brdcst] = linear_separate_window_10thresholds_2H_compbio(wavPS);
    matEcho{t} = echoL;
    matB{t} = brdcst;
    title(sprintf('No. %d threshold',t));
    if t == 1
    [ref1, ref2] = create2Href(brdcst, echoL, Necho, echolength+50); % output E1, E2, the echoes from 1H, and 2H
    elseif g == 1 && t>1
        load('reference12H.mat');
    end
    [de(t).D1, de(t).D2, b] = dechirping(brdcst, echoL, ref1,ref2,wavPS.Necho, echolength, plotflag);
    
end

% Get dechirped image for two harmonics bat calls
% load testmat_bat2H_dechirped2Hs_bw2k



% reshape the "de" struct array into structure that spreads out as echo from
% echo, each sheet is 10 thresholds of one echo
dh1 = struct;
data = zeros(151,10); % one echo, across 10 thresholds  
for j1 = 1:10 % echo
    for h2 = 1:10 % threshold
  
    data(:,h2) = de(h2).D1(:,j1);
    end
    dh1(j1).data = data;
end

dh2 = struct;
data = zeros(151,10);
for k1 = 1:10 % echo
    for s1 = 1:10 % threshold
        data(:,s1) = de(s1).D2(:,k1);
    end
    dh2(k1).data = data;
end
    
[index, pos] = measureWidth(ts, Necho, 500);
% For those with glints longer than 500, change wav parameter and detect
% both of the separated echoes

% 1. find the index
long_in = setdiff(1:Necho, index);
LA = struct;
comp = zeros(length(long_in),1);

% 2. look for the dechirped representations for both
for j = 1:length(long_in)
    tj = long_in(j);
    % form a unique broadcast + echo pair
    comp(j) = pos(tj+1) - 500;
    brd = ts.data(1:pos(1)+ echolength*4);
    echo = ts.data(pos(tj+1) - 500: pos(tj+1)+echolength*4);
    ns = struct;
    ns.data = [zeros(500,1); brd; 1.2*echo];
    ns.fs = ts.fs;

    figure
    SIML2 = runBiscatMain2_noplot(cfg, ns);
    wavPS.simStruct = SIML2;
    wavPS.Necho = 2;
    % find the intercepts
    wavPS.startingThPercent = 10;
    wavPS.whenBrStart = 1E-3;
    wavPS.SepbwBRand1stEchoinSmpls = round(6E-3*ns.fs); % sepearation time between broadcast and 1st echo, in samples
    wavPS.callLenForMostFreq = .4E-3; % twi,
    wavPS.callLenForHighFreq = .4E-3; % twh
    wavPS.callLenSpecial = .4E-3; % twu
    for t1 = 1:wavPS.NT
       
        wavPS.NoT = t1;
        [echoL,brdcst] = linear_separate_window_10thresholds_2H_compbio(wavPS);
        [D1, D2] = dechirp_detached(brdcst, echoL, 2);
         LA(j).d1(:,:,t1) = D1;
         LA(j).d2(:,:,t1) = D2;

    end


end

frange = splitRange(b);
ip_LH1 = findnotches_multiechoes_bat1H(Necho, index, dh1, frange(1,:));
ip_LH2 = findnotches_multiechoes_bat2H(Necho, index, dh2, frange(2,:));




%% Plot 3D SCAT
Fc = 25E3:.5E3:100E3;
des = 1:3:151;
opind = 4*ones(Necho,1);
A = zeros(151,Necho,2); % two harmonics added to the third dimension
for th = 2:4 % plot the 3D-SCAT figure from 2nd threshold to 4th threshold

    for m = 1:Necho
        A(:,m,1) = de(th).D1(:,m)/500; 
        A(:,m,2) = de(th).D2(:,m)/500;
    end
    Nbin = 5*ones(Necho,1);



figure
intvl = .5E3;
intvl_thin = 2E3;
Sep = 2;
fbarea = [1:5 46:58];
CC = flipud(colormap(winter(3*length(des))));
for w = 1:length(index) % echo wise
  
ipforall = [ip_LH1{w} ip_LH2{w}];
IP1 = ip_LH1{w};
kip = ismember(ipforall,fbarea);% eliminate the overlapping area 46:58
k_1H = ismember(IP1,fbarea);
IP1 = IP1(~k_1H);
add1Hscatter(A(:,:,1), w, des, Sep, Fc, CC, .8);
add2Hscatter(A(:,:,2), w, des, Sep, Fc, CC, .3);
threeDplotfunc_cr_bat(A(:,:,1), w, des, intvl,index, IP1, ipforall(~kip), Sep, Fc);
fc = 25:.5:100;
addOverlapPatch(fc(46), fc(58), 25);
addDelayEst(A(:,w,1), A(:,w,2), w, Sep);
end


for q = 1:length(long_in)
    delay = zeros(2,1);
    for i = 1:2
        minD = comp(q)/500+min(LA(q).d1(:,1,th)/500);
        sep = Sep*(q+w-1) - minD;
        add1Hscatter2(comp(q)/500+LA(q).d1(:,i,th)/500, des, sep, Fc, CC, .8);
        add2Hscatter2(comp(q)/500+LA(q).d2(:,i,th)/500, des, sep, Fc, CC, .3);
        delay(i) = addDelayEst_detached(comp(q)/500+LA(q).d1(:,i,th)/500, comp(q)/500+LA(q).d2(:,i,th)/500, sep);
        
    end
    gs = diff(delay);
    str = sprintf('%d \\mus', round(gs*1000));
    text(min(comp(q)/500+LA(q).d1(:,i,th)/500)+sep, 100, 20, str, 'FontSize', 12, 'FontName', 'times');
    
    
end

adjfig_3dscat(th);

end

end