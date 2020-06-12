% The main program to run for the binaural thing
% Created by Chen Ming on 4/2/2019
% this program should be used when the bat is moving and the target is
% stationary
clc;clear;close all;
sandbox;
% load config file
load config_binaural.mat

%%% Fix the bat's head at origin
d = 0.014; % the distance between two ears
coord = [-d/2, 0; d/2, 0]; % initial coordinates of the two ears [xl yl; xr yr]

%%% SOUND SOURCE PARAMETERS -- tars
% load targets structure from 'targets' folder
addpath('binaural-targets');
load('scenario1'); % scenario 1, multiple spacings, randomly generated between 0 - 300 us, except No.8, which is the desired 100 us.

% tars = struct();
% tars.r = [4, 5, 6, 7, 8, 10, 12, 15]; % 8, 10
% tars.theta = [30, 0, -40, -20, -50, 20, 40, -20]; % , 0, 30
% tars.NoG  = 2; % number of glints
% tars.tin = [300, 300, 300, 300, 300, 300, 300, 100]; % , 300, 100
% tars.Ns = length(tars.theta); % number of source
% tars.index = 1:tars.Ns;

ran = tars.index;
scra_range = []; % excluded targets
goal = 100; % this is the target
sense_invl = goal + 100; % initial value, must be a wrong number


gs = []; % glint spacing

rot_step = 15; % rotate 5 degrees one step
alpha = 0; % initial rotation is 0 degrees
vec_ini = [0 1];
flyspeed = 0.5;
coordear_all = zeros(100,4);

vec_all(1,:) = vec_ini;
vec = vec_ini;
tarvec(1) = 1; % the labels of the 
mean_ip = {};

v = 1;
u = 1;

while abs(sense_invl - goal) > 10 % the accuracy of SCAT model

%%% Search From Nearest Target -- Find the nearest target
new_ran = setdiff(ran, scra_range);
ind = findNearestTarget(tars, mean(coord)/2, new_ran);
ind = new_ran(ind);
coord_ss = [tars.r(ind)*cosd(tars.theta(ind)+90), tars.r(ind)*sind(tars.theta(ind)+90)];
delay = [dist(coord_ss, coord(1,:)), dist(coord_ss, coord(2,:))]; % delay in distance [l, r] (m)

%%% Add amplitude - latency trading effect
%off_theta = tars.theta(ind); 

head = mean(coord); % mean: the sum has already been divided by 2
aimvec = coord_ss - head;
off_theta = twovecs(aimvec, vec); % angle off head aim axis, off_theta + : left side to the bat; off_theta - : right side to the bat
oneside_al = 17/60*off_theta*11*1E-6*340; % one side time difference due to amplitude latency trading effect in meter

new_delay = [delay(1)-oneside_al, delay(2)+oneside_al]; % NEW_DELAY is the accurate delay used to calculate the echoes

%%% Calculate the echoes arriving at the bats' two ears.
ts = generate_sigs_with_delay(new_delay);
Fs = ts.fs;
thresh = 20; % threshold to SAMPLE NUMBER difference between two ears
%%% Calculate the time difference received from the two ears by using SCAT
%%% model, initial condition
tsL = struct;
tsL.fs = ts.fs;
tsL.data = ts.data(:,1);

tsR = struct;
tsR.fs = ts.fs;
tsR.data = ts.data(:,2);
figure
SIML = runBiscatMain2_noplot(cfg, tsL);
figure
SIMR = runBiscatMain2_noplot(cfg, tsR);

% set up the wave parameters

wavStructParam = struct;
wavStructParam.Fs = Fs; % sampling frequency of the broadcast-echoes duo
wavStructParam.callLenForMostFreq = .5E-3; % twi,
wavStructParam.callLenForHighFreq = 1E-3; % twh
wavStructParam.callLenSpecial = 1.8E-3; % twu
wavStructParam.sepFlag = 1; % if ==1, use fixed Sep, or determine sep with findSep function when all echoes are so close together, i.e., chain echoes
wavStructParam.whenBrStart = .5E-3; % br_startSec, the time when broadcast started
wavStructParam.startingThPercent = 3; % th_val, integer, the series of thresholds start from 5 percent.
wavStructParam.th_type = 'const'; % th_type
wavStructParam.SepbwBRand1stEchoinSmpls = round(10E-3*Fs); % sepearation time between broadcast and 1st echo, in samples
wavStructParam.color = 'b'; % color setting of the dechirped image.
wavStructParam.ALT = -25; % Coefficient of amplitude latency trading effect, -25 microseconds per dB, must be negative
wavStructParam.NT = 10; % number of thresholds

%%% Get the dechirped response
    figure
for i = 1
    wavStructParam.NoT = i;
    wavStructParam.simStruct = SIML;
    [echoL, firstGapL] = linear_separate_window_10thresholds(wavStructParam);
end
    hold on
    
for j = 1
    wavStructParam.NoT = j;
    wavStructParam.simStruct = SIMR;
    [echoR, firstGapR] = linear_separate_window_10thresholds(wavStructParam);
end

    hold off
figure, h1 = histogram(firstGapL, 'FaceColor', [0.9290 0.6940 0.1250]);
    
[~, max1] = max(h1.Values);
smpl_L = h1.BinEdges(max1);
    
hold on, h2 = histogram(firstGapR, 'FaceColor', [0.3010 0.7450 0.9330]);
legend('Left ear', 'Right ear');
[~, max2] = max(h2.Values);
smpl_R = h2.BinEdges(max2);
    
sense_delay = abs(smpl_R - smpl_L); % in samples

    %%%% Evaluation of flyspeed and rot_step
    if sense_delay > 150 % when target is almost behind the bat
        rot_step = 20;
    elseif sense_delay <80 % when target is basically right in front
        rot_step = 10;
    else
        rot_step = 15;
    end
    
    if smpl_R < round(1.5*2/340*Fs) % 1.5 m
        flyspeed = .2;
    elseif smpl_R > round(3*2/340*Fs) % more than 3m
        flyspeed = .5;
    else
        flyspeed = .4;
    end


%%% Rotate the bat's head/two ears slowly to adjust the position until the
%%% time difference was not distinguishable in the histograms
%%% while loop




while abs(sense_delay) > thresh 
    %v
    coordear_all(v,:) = [coord(1,:) coord(2,:)];
    vec_all(v,:) = vec;
    
    distdiff = diff(new_delay);
    alpha = alpha + sign(distdiff)*rot_step; % counterclockwise: +, clockwise: -
    rot_step_wsign = sign(distdiff)*rot_step;
    R5 = [cosd(rot_step_wsign) -sind(rot_step_wsign); sind(rot_step_wsign) cosd(rot_step_wsign)];
    % rotation matrix for the vector
    R = [cosd(alpha) -sind(alpha); sind(alpha) cosd(alpha)];
    
    vec = vec_ini*R';
    
    % move the bat's head toward the direction of the vector
%     if u <= 2
%         coord = [d/2*cosd(alpha+180), d/2*sind(alpha+180);d/2*cosd(alpha), d/2*sind(alpha)] + flyspeed*vec;
%     else
        % move the coord back to the origin
        xm = sum(coord(:,1))/2;
        ym = sum(coord(:,2))/2;
        coord_x = coord(:,1) - xm;
        coord_y = coord(:,2) - ym;
        % rotate
        coord_rot = [coord_x coord_y]*R5';
        % then move back
        coord_trans = [coord_rot(:,1)+xm coord_rot(:,2)+ym];
        
        coord = coord_trans + flyspeed*vec;

    head = head + flyspeed*vec;
    if size(coord_ss,1) ~= size(head,1)
        coord_ss = coord_ss';
    end
    aimvec = coord_ss - head;
    off_theta = twovecs(aimvec, vec);
    %off_theta = off_theta - alpha; % this comes with a sign
    delay = [dist(coord_ss, coord(1,:)), dist(coord_ss, coord(2,:))];
    if delay(1)>delay(2)
        off_theta = (-1)*abs(off_theta);
    else
        off_theta = abs(off_theta);
    end
    
    oneside_al = 17/60*off_theta*11*1E-6*340;
    
    
    new_delay = [delay(1)-oneside_al, delay(2)+oneside_al];
    
    ts = generate_sigs_with_delay(new_delay);
    %%% Calculate the time difference received from the two ears by using SCAT
    %%% model, initial condition
    tsL = struct;
    tsL.fs = ts.fs;
    tsL.data = ts.data(:,1);

    tsR = struct;
    tsR.fs = ts.fs;
    tsR.data = ts.data(:,2);
    SIML = runBiscatMain2_noplot(cfg, tsL);
    SIMR = runBiscatMain2_noplot(cfg, tsR);
    
    
    for i = 1
        figure
        wavStructParam.NoT = i;
        wavStructParam.simStruct = SIML;
        [echoL, firstGapL] = linear_separate_window_10thresholds(wavStructParam);
    end

    for j = 1
        hold on
        wavStructParam.NoT = j;
        wavStructParam.simStruct = SIMR;
        [echoR, firstGapR] = linear_separate_window_10thresholds(wavStructParam);
        
    end

    figure, h1 = histogram(firstGapL, 'FaceColor', [0.9290 0.6940 0.1250]);
    
    [~, max1] = max(h1.Values);
    smpl_L = h1.BinEdges(max1);
    
    hold on, h2 = histogram(firstGapR, 'FaceColor', [0.3010 0.7450 0.9330]);
    
    [~, max2] = max(h2.Values);
    smpl_R = h2.BinEdges(max2);
    legend('Left ear', 'Right ear');
    sense_delay = abs(smpl_R - smpl_L);
    %coordear_all(u,:) = [coord(1,:) coord(2,:)];
    fprintf('sense delay %d \n', sense_delay);

    close all;
    
    %%%% Evaluation of flyspeed and rot_step
    if sense_delay > 150 % when target is almost behind the bat
        rot_step = 20;
    elseif sense_delay <80 % when target is basically right in front
        rot_step = 10;
    else
        rot_step = 15;
    end
    
    if smpl_R < round(1.5*2/340*Fs) % 1.5 m
        flyspeed = .2;
    elseif smpl_R > round(3*2/340*Fs) % more than 3m
        flyspeed = .5;
    else
        flyspeed = .4;
    end
    
    v = v+1;
    tarvec(v) = ind;
end
    
    %%%%% define new range depends on if this target is the desired one
    % get echo with the characteristics of THE target : ind

    TS = generate_sigs_with_delay_multiglints(dist(head,coord_ss), tars.tin(ind));
    SIMH = runBiscatMain2_noplot(cfg, TS);
    d1 = zeros(81, 10);   
    Fc = SIMH.coch.Fc;
    Nth = 10;
    Mat = cell(Nth, length(Fc));
    close all;
    for t = 1:10
        wavStructParam.NoT = t;
        wavStructParam.simStruct = SIMH;
        [echoR, firstGapR] = linear_separate_window_10thresholds(wavStructParam);
        %[echoH, firstGapH] = linear_separate_window_10thresholds(SIMH, t, Fs, 0.5E-3, 1E-3, 1.8E-3, 15E-3, 5000, 'const', 5, 'b');
        d1(:,t) = firstGapH;
        Mat(t,:) = echoH;
    end

    [row, col] = find(isnan(d1));
    
    mean_ip{u} = findnotches2(d1, col(1));
   
    hp = histogram(diff(mean_ip{u}));
    [~, maxp] = max(hp.Values);
    frq_intl = (hp.BinEdges(maxp)+hp.BinEdges(maxp+1))/2;
    sense_invl = 1/frq_intl*1E3;
    gs = [gs sense_invl]; % glint spacing estimates for the targets
    u = u + 1;
    %%%% if sense_intl is not desired, then scratch 
    if abs(sense_invl - goal) > 10
        scra_range = [scra_range ind];
    end
    fprintf('on target %d \n', ind);
    fprintf('the glint spacing is %d, the difference is %d, true glint spacing %d \n', round(sense_invl), round(sense_invl)-100, tars.tin(ind));
end
    vec_all(v,:) = vec; % the last one
    coordear_all(v,:) = [coord(1,:) coord(2,:)];
    %%% fly the bat to the target after turning/ or receiving the echo from the
    %%% target after the bat is facing the right target
% 1. get the echo with two glints for example with time interval 100us
% 2. use scat model to sense the delay
% 3. Decide if this is the right target


%%% Approach with the right direction
%%% if steps are too small, approach without rotation
%%% The bat should move toward the right direction before it reaches the
%%% target, it's the capability

% function [xx, yy, zz] = generateTarget(numOfTarget)
% % To generate certain number of targets at different locations
% % input, numOfTarget, number of targets
% 
% 
% end

function I = findNearestTarget(tarStruc, batpos, range)
% Find the nearest target by inputing the target structure: tarStruc
% bat position, batpos, which is the head position
% range, is the indices of the targets that are not examined yet.
    tr = tarStruc.r(range);
    tthe = tarStruc.theta(range) + 90;
    tcoord_x = tr.*cosd(tthe);
    tcoord_y = tr.*sind(tthe);
    tc = [tcoord_x'  tcoord_y'];
    alldist = distV(tc, batpos);
    
    [~, I] = min(alldist);
end

function dd = dist(a,b)
    dd = sqrt((a(1) - b(1))^2 + (a(2) - b(2))^2);
end

function dis = distV(a,b)
    % a is a vector [x y], b is a point
    x = a(:,1) - b(1);
    y = a(:,2) - b(2);
    
    dis = sqrt(x.^2 + y.^2);
end

function ang = twovecs(u,v)
% angle in between two vectors in degrees
    u = [u 0];
    v = [v 0];
    ang = atan2(norm(cross(u,v)),dot(u,v));
    ang = rad2deg(ang);
    r = vrrotvec(v,[0,1,0]);
    m = vrrotvec2mat(r);
    if size(u,1)<size(u,2)
        newu = u*m';
    else
        newu = u*m;
    end
    if newu(1)>0
        ang = -ang;
    end
%     w = u - v;
%     if w(1)>0
%         ang = -ang;
%     end
end

