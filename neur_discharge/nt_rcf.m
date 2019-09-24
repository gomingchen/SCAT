function neur = nt_rcf(fs, bmm, Fc, cfg)
% NT_RCF  Runs Rectify, Compress, and Filter neural model
%
% neur = NT_RCF(FS,BMM,FC,CFG)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rectification
switch cfg.neur_rcf_rect
    case 1              % half-wave
        rcf = bmm;
        rcf(rcf<0)=0;
    case 2              % full-wave
        rcf = abs(bmm);
    case 3              % bypass
        rcf = bmm;
    otherwise
        warning('NEURTRAN:RCF','unknown rectification parameter')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compression
switch cfg.neur_rcf_comp
    %  Need to apply cochlear voltage response curve (Fig. 3.6, pg. 50) to
    %    include non-linearities (logarithmic is a good assumption).
    %  This is necessary to boost low amplitude onset responses so Neural
    %    Discharge stage is more accurate for time estimates.
    case 1              % logarithmic
        rcf = log10(rcf);
        rcf(rcf == -Inf) = -120;
    case 2              % bypass

    otherwise
        warning('NEURTRAN:RCF','unknown compression parameter')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering (Envelope Smoothing)
switch cfg.neur_rcf_filt
    case 1              % Butterworth LPF
        [neur.filt_b, neur.filt_a] = butter(cfg.neur_rcf_ord, cfg.neur_rcf_fc/(fs/2));
        rcf = filter(neur.filt_b, neur.filt_a, rcf, [], 1);
    case 2              % bypass
        neur.filt_b = 1;
        neur.filt_a = 1;
    otherwise
        warning('NEURTRAN:RCF','unknown filtering parameter')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike Generation

% find refractory periods in samples
Na = round(1e-3*cfg.neur_rcf_refaper*fs);
Nr = round(1e-3*cfg.neur_rcf_refrper*fs);

% construct anonymous function for the specified relative refractory period
switch cfg.neur_rcf_reffcn
    case 1
        % constant threshold
        refFun = @(val,loc) val(1) * ones(1,length(loc));
    case 2
        % linearly decreasing threshold
        refFun = @(val,loc) -val(1)/Nr * loc + val(1);
    case 3
        % exponentially decreasing threshold
        refFun = @(val,loc) val(1) * exp(-2.3*loc/Nr);
    otherwise
        error('Invalid selection for refractory period')
end

% detect events of interest and fire spikes
sVal = cell(cfg.coch_steps,1);
sLoc = cell(cfg.coch_steps,1);
switch cfg.neur_rcf_spike
    % Fixed threshold
    case 1
        thresh = sqrt(mean(rcf.^2,1));   % just use RMS value of each data series for now
        for i=1:cfg.coch_steps
            [sVal{i}, sLoc{i}] = findpeaks(rcf(:,i), 'MINPEAKHEIGHT', thresh(i));
        end

    % Variable threshold
    case 2
        % Threshold crossing afferent neurons mark an event whenever the
        % level of the signal crosses a predetermined threshold for a given
        % filter band.  The firing thresholds are determined adaptively
        % using the known broadcast SNR

        thresh = zeros(1,cfg.coch_steps);
        for i=1:cfg.coch_steps
            [sVal{i}, sLoc{i}] = findpeaks(rcf(:,i), 'THRESHOLD', thresh(i));
        end
        % find (diff(sVal{i}) > thresh/time)

    % Peak detection neurons
    case 3
        % Peak-detection afferent neurons fire only when a peak is above a
        % predetermined threshold.  During the refractory period the firing
        % threshold is increased to the current peak's level.

        for i=1:cfg.coch_steps
            % find all maximums in time series within +/- 1/4 period from largest peaks
            [sVal{i}, sLoc{i}] = findpeaks(rcf(:,i),'MINPEAKDISTANCE',round(fs/(Fc(i)*4)));

            % remove all spikes not meeting minimum threshold
            idx = find(sVal{i} > cfg.neur_rcf_thrlev);
            [sVal{i}, sLoc{i}] = deal(sVal{i}(idx), sLoc{i}(idx));
            
            % thin spikes using refractory function
            [sVal{i}, sLoc{i}] = thinspikes(sVal{i}, sLoc{i}, refFun, Nr);
        end

    % Bypass spike generation (useful for matched filt, xcorr models, etc.)
    case 4

    otherwise
        warning('NEURTRAN:RCF','unknown spike threshold parameter')
        return
end

% return RCF output as "probability of spike"
neur.prob = rcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% record spike trains to output struct
neur.spikes = sparse(size(rcf,1),size(rcf,2),false);  %initialize sparse neural spike output
for i=1:cfg.coch_steps
    neur.spikes(sLoc{i},i) = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rVal,rLoc,thresh] = thinspikes(sVal,sLoc,refFun,N)
% Thin spikes according to function handle rfFun

% init threshold vector using refractory function
thresh = zeros(length(sLoc),1);

% init vector of thinned spikes to be returned
rLoc = zeros(length(sLoc),1);
rVal = zeros(length(sLoc),1);
rIdx = 1;        % index of return vectors

% raise threshold by refractory period for each consecutive valid peak
for j=1:length(sVal);
    
    % test current spike against threshold
    if sVal(j) > thresh(j)
        
        % mark spike index as valid
        rLoc(rIdx) = sLoc(j);
        rVal(rIdx) = sVal(j);
        rIdx = rIdx+1;
        
        % update threshold function, if necessary
        k = find(sLoc < sLoc(j)+N, 1, 'last');
        thresh(j:k) = refFun(sVal(j:k),sLoc(j:k));
        
    end
end

% truncate unused memory
rLoc(rIdx:end) = [];
rVal(rIdx:end) = [];
